# -*- coding: utf-8 -*-
"""
TesseraHabitat.pyt

Habitatprediktion ur Tessera-embeddings: utgå från kända fyndpunkter för en art,
hämta embedding-vektorn i varje fyndpunkt, och måla ut ett likhetsraster över ett
sökområde så att platser som liknar de kända lokalerna framträder.

Vilka tiles som finns, var de ligger och vilka år som är publicerade kommer från
biblioteket geotessera, som också hämtar landmaskerna. Själva fyndpunkterna läses
däremot med Range-anrop direkt ur Tesseras filer, inte via geotessera: en tile är
ca 90 MB och geotessera hämtar alltid hela tiles, medan fynddata för en art ligger
utspridda över en hel region. Ett fynd kostar två små anrop i stället för 90 MB,
och den skillnaden är hela förutsättningen för verktyget.

Filformatet är okomprimerad .npy, så en enskild pixel kan läsas ut med ett
byte-intervall:

    .../{år}/grid_{lon}_{lat}/grid_{lon}_{lat}.npy
        .npy-huvud (128 byte) + int8 i C-ordning, form (H, W, 128)
        pixelns vektor = 128 byte på offset huvud + (rad * W + kolumn) * 128

    .../{år}/grid_{lon}_{lat}/grid_{lon}_{lat}_scales.npy
        .npy-huvud (128 byte) + float32 i C-ordning, form (H, W)
        pixelns skalfaktor = 4 byte på offset huvud + (rad * W + kolumn) * 4

    landmasken för tilen, 13,6 kB, hämtas via geotessera och är den enda
    källan till tilens koordinatsystem och origo

Embeddingen är kvantiserad som int8 med en skalfaktor per pixel — en enda skalär
för alla 128 kanaler — så det verkliga värdet är int8 * scale. Kostnaden per
fyndpunkt blir två små anrop, plus en landmask per tile som återanvänds ur cachen.

Uppmätt egenhet hos servern som koden måste ta hänsyn till:

  * Flerdelade Range-huvuden ("bytes=a-b, c-d") ignoreras: servern svarar 200 med
    hela filen i stället för 206 med delarna. Varje läsning kontrollerar därför
    att statuskoden är 206 innan kroppen läses, och varje sammanhängande löpa
    hämtas med ett eget anrop. Utan den kontrollen blir en avsedd 128-byteläsning
    en nedladdning av hela tilen. Kontrollerat mot S3-adressen; kravet på
    User-Agent gällde den tidigare spegeln på data.source.coop och finns inte här,
    men huvudet skickas ändå.

Likhetsberäkningen är densamma som i "Tessera similarity search" i denna mapp:
skalärprodukten mellan pixelns vektor och referensvektorn, med normering ger det
kosinuslikhet mellan -1 och 1.

    https://developers.google.com/earth-engine/tutorials/community/satellite-embedding-05-similarity-search

Skillnaden mot det verktyget är hur många referenser som vägs samman och hur.
En art använder ofta flera habitat, och fynddata är därmed multimodala i
embedding-rummet. Standardvalet här är därför "största likhet mot någon
fyndpunkt": varje pixel jämförs med samtliga fynd och behåller sin bästa
träff, vilket bevarar flera habitat i stället för att medelvärdesbilda bort dem.

Källa : https://geotessera.org/  (data via s3://tessera-embeddings)
Krav  : ArcGIS Pro 3.x (arcpy), numpy och geotessera.
        geotessera finns inte i standardmiljön arcgispro-py3. Klona miljön och
        installera paketet, t.ex. arcgispro-py3-personal, och peka Pro på den
        kloningen innan verktygslådan används.
"""

import ast
import concurrent.futures
import math
import os
import re
import shutil
import tempfile
import time
import urllib.error
import urllib.request

import numpy as np

import arcpy

try:
    from geotessera import GeoTessera
    from geotessera import registry as gt_registry
    _GEOTESSERA_ERROR = None
except Exception as _exc:                                   # noqa: BLE001
    GeoTessera = None
    gt_registry = None
    _GEOTESSERA_ERROR = _exc

# ── Konstanter ────────────────────────────────────────────────────────────────

WGS84_WKID = 4326

TILE_DEG = 0.1          # tile-sida i grader
N_CHANNELS = 128        # kanaler per pixel i Tessera-embeddingen

# Dataset-versioner: etikett -> (version, variant) som geotessera vill ha dem.
# Samma tabell som i "Tessera embeddings to GDB".
DATASETS = {
    "v1 (global, 2017-2025)":      ("v1", "vultr"),
    "v1.1 Cambridge (regional)":   ("v1.1", "cambridge"),
    "v2 beta (delvis täckning)":   ("v2", "2B-L~beta1"),
}
DEFAULT_DATASET = "v1 (global, 2017-2025)"

YEARS = [str(y) for y in range(2017, 2026)]
DEFAULT_YEAR = "2024"

AGG_MAX = "Största likhet mot någon fyndpunkt (behåller flera habitat)"
AGG_MEAN = "Medelvektor av alla fyndpunkter (ett habitat)"

CAT_TESSERA = "Tessera-data"
CAT_REF = "Referens"
CAT_OUT = "Utdata"
CAT_MAP = "Karta"

STATUS_USED = "använd"
STATUS_DROPPED = "utesluten"
STATUS_NODATA = "saknar data"
STATUS_NOTILE = "ingen tile"

_USER_AGENT = "TesseraHabitat.pyt (ArcGIS Pro)"
_HTTP_TIMEOUT = 120
_HTTP_ATTEMPTS = 4
_RETRY_STATUS = (408, 429, 500, 502, 503, 504)
_READ_CHUNK = 1024 * 1024
_WORKERS = 8

_CACHE_DIRNAME = "Tessera_nedladdning"
_SCRATCH_DIRNAME = "TesseraHabitat_arbetsmapp"

# Mappar som synkas till molnet — olämpliga som cache.
_SYNC_HINTS = ("onedrive", "sharepoint", "dropbox", "google drive")

# Minnesbudget per inläst radblock i likhetsberäkningen (bytes). Både
# indataläsningen och matrisen med likheter mot samtliga referensvektorer
# ryms inom den.
_BLOCK_BUDGET_BYTES = 256 * 1024 ** 2

# Minnesbudget per block i grannskapsberäkningen mellan fyndpunkter.
_KNN_BUDGET_BYTES = 64 * 1024 ** 2

# Varna om utdatarastret blir större än så här (bredd x höjd), eftersom hela
# likhetsrastret hålls i minnet som en float32-array.
_MAX_SANE_PIXELS = 150_000_000

# Största tillåtna samplingsradie runt en fyndpunkt (m).
_MAX_RADIUS_M = 500.0


class _TileMissing(Exception):
    """Tilen finns inte publicerad (HTTP 404). Vanligt över öppet vatten."""


def _sr(wkid):
    return arcpy.SpatialReference(wkid)


def _sr_is_valid(sr):
    """
    Är sr ett användbart koordinatsystem?

    factoryCode duger inte ensamt som test: ett eget definierat koordinatsystem
    har koden 0 men är fullt giltigt, medan ett tomt SpatialReference också har
    koden 0. Det som skiljer dem är att det tomma saknar WKT-definition.
    """
    if sr is None:
        return False
    try:
        if sr.factoryCode:
            return True
        return bool(sr.exportToString())
    except Exception:
        return False


def _project_geometry(geom, target_sr):
    """
    Projicera en verklig geometri (ur en featureklass) till target_sr. En sådan
    geometri bär alltid sitt eget koordinatsystem, till skillnad från ett
    GPExtent-objekt, så projectAs kan användas direkt.
    """
    try:
        return geom.projectAs(target_sr)
    except Exception as exc:
        raise ValueError(
            "Kunde inte omvandla en fyndpunkt till {}: {}".format(target_sr.name, exc)
        )


# =============================================================================
# Tile-rutnät och adresser
# =============================================================================

def _tile_center(index):
    """Tile-centrum för ett heltalsindex: index 180 -> 18.05."""
    return round(index * TILE_DEG + TILE_DEG / 2.0, 2)


def _tile_of(lon, lat):
    """Tile-centrum (lon, lat) för punkten (lon, lat) i WGS84."""
    return (_tile_center(int(math.floor(lon * 10))),
            _tile_center(int(math.floor(lat * 10))))


def _grid_name(lon, lat):
    """Filnamnsstammen för en tile, t.ex. 'grid_18.05_59.35'."""
    return "grid_{:.2f}_{:.2f}".format(lon, lat)


def _require_geotessera():
    """Ge ett begripligt fel när paketet saknas i den aktiva Python-miljön."""
    if GeoTessera is None:
        import sys
        raise ValueError(
            "Paketet geotessera kunde inte laddas i den Python-miljö som ArcGIS Pro "
            "använder ({}). Klona arcgispro-py3, installera geotessera i kloningen "
            "och byt aktiv miljö i Pro under Settings, Package Manager. "
            "Ursprungligt fel: {}".format(
                os.path.basename(os.path.normpath(sys.prefix)), _GEOTESSERA_ERROR)
        )


# En klient per (version, variant, cache-mapp). Att öppna registret läser ett
# manifest över samtliga publicerade tiles, så klienten återanvänds i sessionen.
_client_cache = {}


def _client(dataset, cache_dir, messages=None):
    """GeoTessera-klient för en dataset-etikett, med landmasker i cache_dir."""
    _require_geotessera()
    if dataset not in DATASETS:
        raise ValueError("Okänd dataset-version: {}".format(dataset))
    version, variant = DATASETS[dataset]

    key = (version, variant, os.path.abspath(str(cache_dir)))
    if key in _client_cache:
        return _client_cache[key]

    if messages is not None:
        messages.addMessage(
            "Läser Tessera-registret ({} {})... första gången hämtas ett "
            "manifest över alla tiles, vilket tar en stund.".format(version, variant)
        )
    try:
        client = GeoTessera(dataset_version=version, dataset_variant=variant,
                            embeddings_dir=str(cache_dir))
    except Exception as exc:                                # noqa: BLE001
        raise ValueError(
            "Kunde inte läsa Tessera-registret för {} {}: {}".format(version, variant, exc)
        )
    _client_cache[key] = client
    return client


def _dataset_base(client, year):
    """
    Bas-URL för en tiles filer, byggd ur geotesseras egna konstanter.

    Adressen härleds hellre än hårdkodas: geotessera 0.9 flyttade datat från
    data.source.coop till S3, och den flytten ska slå igenom här utan att
    sökvägar skrivs av för hand.
    """
    root = gt_registry.TESSERA_BASE_URL.rstrip("/")
    version = client.registry._version_path if hasattr(client.registry, "_version_path")         else DATASETS[DEFAULT_DATASET][0]
    return "{}/{}/{}/{}".format(root, version, gt_registry.EMBEDDINGS_DIR_NAME, year)


def _embedding_url(client, year, lon, lat):
    name = _grid_name(lon, lat)
    return "{}/{}/{}.npy".format(_dataset_base(client, year), name, name)


def _scales_url(client, year, lon, lat):
    name = _grid_name(lon, lat)
    return "{}/{}/{}_scales.npy".format(_dataset_base(client, year), name, name)


# =============================================================================
# HTTP
# =============================================================================

def _request(url, start=None, length=None):
    """
    En begäran med User-Agent, och med Range-huvud när ett intervall begärs.

    Bara ett enda intervall per begäran: servern struntar i flerdelade
    Range-huvuden och skickar hela filen i stället, kontrollerat mot S3.
    User-Agent krävdes av den tidigare spegeln på data.source.coop, som svarade
    403 utan huvudet. S3 bryr sig inte, men huvudet skickas ändå så att anropen
    går att känna igen i serverloggar.
    """
    headers = {"User-Agent": _USER_AGENT}
    if start is not None:
        headers["Range"] = "bytes={}-{}".format(start, start + length - 1)
    return urllib.request.Request(url, headers=headers)


def _read_exactly(resp, length):
    parts = []
    got = 0
    while got < length:
        chunk = resp.read(min(_READ_CHUNK, length - got))
        if not chunk:
            break
        parts.append(chunk)
        got += len(chunk)
    return b"".join(parts)


def _http_bytes(url, start=None, length=None):
    """
    Hämta hela filen, eller exakt length byte från offset start.

    Tillfälliga fel görs om med exponentiell backoff. TimeoutError är en
    OSError och HTTPError ärver både URLError och OSError, så except-satserna
    måste stå i den ordning de gör här.
    """
    for attempt in range(_HTTP_ATTEMPTS):
        last = attempt == _HTTP_ATTEMPTS - 1
        try:
            with urllib.request.urlopen(
                    _request(url, start, length), timeout=_HTTP_TIMEOUT) as resp:
                if start is not None and getattr(resp, "status", None) != 206:
                    # Servern ignorerade Range och skulle skicka hela filen —
                    # kroppen läses medvetet inte, den kan vara 90 MB.
                    raise IOError(
                        "Servern besvarade inte Range-begäran mot {} med 206.".format(url)
                    )
                data = _read_exactly(resp, length) if length is not None else resp.read()

            if length is not None and len(data) != length:
                raise IOError(
                    "Fick {} av {} byte från {}.".format(len(data), length, url)
                )
            return data

        except urllib.error.HTTPError as exc:
            if exc.code == 404:
                raise _TileMissing(url)
            if exc.code not in _RETRY_STATUS or last:
                raise
        except (urllib.error.URLError, OSError):
            if last:
                raise
        time.sleep(2 ** attempt)

    raise IOError("Kunde inte hämta {}.".format(url))


def _http_error_msg(exc):
    if isinstance(exc, urllib.error.HTTPError):
        return "HTTP-fel {} från Tessera-servern: {}".format(exc.code, exc.reason)
    if isinstance(exc, urllib.error.URLError):
        return ("Kunde inte nå Tessera-servern ({}). Kontrollera internetanslutning "
                "och eventuell proxy.".format(exc.reason))
    if isinstance(exc, TimeoutError):
        return "Servern svarade inte i tid. Försök igen."
    return "Fel vid anrop till Tessera-servern: {}".format(exc)


def _npy_header(url):
    """
    (dataoffset, form) ur ett .npy-huvud.

    Huvudet är 128 byte i de filer Tessera publicerar i dag, men längden läses
    ur filen i stället för att antas — det kostar inget extra anrop eftersom
    de första 256 byten hämtas i ett svep ändå.
    """
    raw = _http_bytes(url, 0, 256)
    if raw[:6] != b"\x93NUMPY":
        raise ValueError("{} är inte en .npy-fil.".format(url))
    if raw[6] == 1:
        header_len = int.from_bytes(raw[8:10], "little")
        offset = 10
    else:
        header_len = int.from_bytes(raw[8:12], "little")
        offset = 12
    if offset + header_len > len(raw):
        raw = _http_bytes(url, 0, offset + header_len)
    info = ast.literal_eval(raw[offset:offset + header_len].decode("latin1"))
    if info.get("fortran_order"):
        raise ValueError("{} är lagrad i Fortran-ordning och kan inte läsas styckvis.".format(url))
    return offset + header_len, tuple(info["shape"])


def _fetch_landmask(client, lon, lat):
    """
    Landmasken för en tile. geotessera hämtar den en gång till cache-mappen och
    återanvänder den sedan, så den är billig att be om per fyndpunkt.
    """
    try:
        return client.registry.fetch_landmask(lon=lon, lat=lat)
    except Exception as exc:                                # noqa: BLE001
        raise _TileMissing(
            "Landmasken för {} kunde inte hämtas: {}".format(_grid_name(lon, lat), exc)
        )


# =============================================================================
# Punktsampling ur Tessera
# =============================================================================

class _Observation:
    """En fyndpunkt på väg genom körningen."""

    __slots__ = ("index", "oid", "shape", "tile", "row", "col",
                 "vector", "status", "neighbour", "score", "clipped")

    def __init__(self, index, oid, shape, tile):
        self.index = index
        self.oid = oid
        self.shape = shape
        self.tile = tile
        self.row = None
        self.col = None
        self.vector = None
        self.status = STATUS_NODATA
        self.neighbour = None
        self.score = None
        self.clipped = False


class _TileSource:
    """
    Adresserna och måtten som krävs för att läsa enskilda pixlar ur en tile.

    Georefereringen kommer ur landmasken: .npy-filerna innehåller bara
    pixelvärden och vet ingenting om var de ligger.
    """

    __slots__ = ("lon", "lat", "emb_url", "sca_url", "emb_offset", "sca_offset",
                 "height", "width", "sr", "extent", "cell_w", "cell_h", "landmask")

    def __init__(self, client, year, lon, lat):
        self.lon = lon
        self.lat = lat
        self.emb_url = _embedding_url(client, year, lon, lat)
        self.sca_url = _scales_url(client, year, lon, lat)

    def load_remote(self, client):
        """Nätverksdelen: landmask till cachen och båda .npy-huvudena.
        Innehåller inga arcpy-anrop och kan därför köras i en trådpool."""
        self.landmask = _fetch_landmask(client, self.lon, self.lat)
        self.emb_offset, emb_shape = _npy_header(self.emb_url)
        self.sca_offset, sca_shape = _npy_header(self.sca_url)

        if len(emb_shape) != 3 or emb_shape[2] != N_CHANNELS:
            raise ValueError(
                "{} har formen {}, förväntade (H, W, {}).".format(
                    _grid_name(self.lon, self.lat), emb_shape, N_CHANNELS)
            )
        if tuple(sca_shape[:2]) != tuple(emb_shape[:2]):
            raise ValueError(
                "Skalfilen för {} matchar inte embeddingen.".format(
                    _grid_name(self.lon, self.lat))
            )
        self.height, self.width = emb_shape[0], emb_shape[1]

    def load_georeference(self):
        """arcpy-delen: koordinatsystem, origo och cellstorlek ur landmasken.
        Körs på huvudtråden, arcpy är inte trådsäkert."""
        raster = arcpy.Raster(self.landmask)
        self.sr = raster.spatialReference
        self.extent = raster.extent
        self.cell_w = raster.meanCellWidth
        self.cell_h = raster.meanCellHeight
        if not _sr_is_valid(self.sr):
            raise ValueError(
                "Landmasken för {} saknar koordinatsystem.".format(
                    _grid_name(self.lon, self.lat))
            )
        if (raster.height, raster.width) != (self.height, self.width):
            raise ValueError(
                "Tilen {} matchar inte sin landmask ({}x{} mot {}x{}).".format(
                    _grid_name(self.lon, self.lat), self.width, self.height,
                    raster.width, raster.height)
            )

    def rowcol(self, x, y):
        """Rad och kolumn för en koordinat i tilens eget koordinatsystem."""
        col = int((x - self.extent.XMin) / self.cell_w)
        row = int((self.extent.YMax - y) / self.cell_h)
        return row, col

    def read_run(self, row, col, count):
        """
        En sammanhängande löpa om count pixlar på rad row, från kolumn col.

        Två anrop: ett för int8-vektorerna (count * 128 byte i följd, eftersom
        arrayen ligger i C-ordning) och ett för skalfaktorerna.
        """
        flat = row * self.width + col
        quantized = np.frombuffer(
            _http_bytes(self.emb_url,
                        self.emb_offset + flat * N_CHANNELS,
                        count * N_CHANNELS),
            dtype=np.int8,
        ).reshape(count, N_CHANNELS)
        scales = np.frombuffer(
            _http_bytes(self.sca_url, self.sca_offset + flat * 4, count * 4),
            dtype="<f4",
        )
        return quantized, scales


def _window(source, x, y, radius):
    """
    Cellerna som ska samplas för en punkt: bara mittcellen när radien är noll,
    annars alla celler vars centrum ligger inom radien.

    Returnerar (rader, per_rad, klippt) där rader är en lista med
    (rad, första_kolumn, antal) och per_rad en lista med boolska masker.
    Ett fönster som sticker ut över tilekanten klipps mot tilen — den
    angränsande tilens del hämtas inte, och att det skett rapporteras.
    """
    row, col = source.rowcol(x, y)
    if not (0 <= row < source.height and 0 <= col < source.width):
        return None, None, False

    if radius <= 0:
        return [(row, col, 1)], [np.ones(1, dtype=bool)], False

    span = int(math.ceil(radius / min(source.cell_w, source.cell_h)))
    r0, r1 = row - span, row + span + 1
    c0, c1 = col - span, col + span + 1
    clipped = r0 < 0 or c0 < 0 or r1 > source.height or c1 > source.width
    r0, r1 = max(0, r0), min(source.height, r1)
    c0, c1 = max(0, c0), min(source.width, c1)

    runs, masks = [], []
    xs = source.extent.XMin + (np.arange(c0, c1) + 0.5) * source.cell_w
    for r in range(r0, r1):
        cy = source.extent.YMax - (r + 0.5) * source.cell_h
        inside = ((xs - x) ** 2 + (cy - y) ** 2) <= radius ** 2
        if inside.any():
            runs.append((r, c0, c1 - c0))
            masks.append(inside)
    if not runs:
        return [(row, col, 1)], [np.ones(1, dtype=bool)], clipped
    return runs, masks, clipped


def _sample_one(source, obs, x, y, radius):
    """Embedding-vektorn för en fyndpunkt, dekvantiserad och medelvärdesbildad
    över fönstret. Returnerar None när ingen cell hade giltig data."""
    runs, masks, clipped = _window(source, x, y, radius)
    obs.clipped = clipped
    if runs is None:
        return None

    total = np.zeros(N_CHANNELS, dtype=np.float64)
    n_valid = 0
    for (row, col, count), mask in zip(runs, masks):
        quantized, scales = source.read_run(row, col, count)
        # Skalfaktorn är en enda skalär per pixel, gemensam för alla 128
        # kanaler. NoData markeras som +inf i skalfilen.
        valid = mask & np.isfinite(scales) & (scales > 0)
        if not valid.any():
            continue
        total += (quantized[valid].astype(np.float64)
                  * scales[valid][:, np.newaxis]).sum(axis=0)
        n_valid += int(valid.sum())

    if not n_valid:
        return None
    return (total / n_valid).astype(np.float32)


def _read_observations(points, messages):
    """Fyndpunkterna som _Observation, med tile-tillhörighet i WGS84."""
    wgs84 = _sr(WGS84_WKID)
    described = arcpy.Describe(points)
    if described.shapeType != "Point":
        raise ValueError(
            "Fyndlagret måste innehålla punkter (har {}).".format(described.shapeType)
        )
    if not _sr_is_valid(described.spatialReference):
        raise ValueError("Fyndlagret saknar koordinatsystem.")

    observations = []
    empty = 0
    with arcpy.da.SearchCursor(points, ["OID@", "SHAPE@"]) as cursor:
        for index, (oid, shape) in enumerate(cursor):
            if shape is None:
                empty += 1
                continue
            geographic = _project_geometry(shape, wgs84)
            point = geographic.centroid
            observations.append(
                _Observation(index, oid, shape, _tile_of(point.X, point.Y))
            )

    if empty:
        messages.addWarningMessage(
            "{} fyndpunkt(er) saknade geometri och hoppades över.".format(empty)
        )
    if not observations:
        raise ValueError("Fyndlagret innehåller inga punkter.")
    return observations


def _sample_observations(observations, client, year, radius, messages):
    """
    Hämta embedding-vektorn för varje fyndpunkt med Range-anrop.

    Körningen sker i tre steg eftersom arcpy inte är trådsäkert: allt som bara
    är HTTP görs parallellt, allt som rör arcpy görs på huvudtråden.
    """
    tiles = sorted({obs.tile for obs in observations})
    messages.addMessage(
        "{} fyndpunkt(er) i {} tile(s). Hämtar georeferering...".format(
            len(observations), len(tiles))
    )

    # Steg 1, parallellt: landmask till cachen och .npy-huvuden per tile.
    sources = {}
    missing = set()
    with concurrent.futures.ThreadPoolExecutor(max_workers=_WORKERS) as pool:
        futures = {}
        for lon, lat in tiles:
            source = _TileSource(client, year, lon, lat)
            futures[pool.submit(source.load_remote, client)] = (lon, lat, source)
        for future in concurrent.futures.as_completed(futures):
            lon, lat, source = futures[future]
            try:
                future.result()
            except _TileMissing:
                missing.add((lon, lat))
                continue
            except (urllib.error.URLError, OSError) as exc:
                raise ValueError(_http_error_msg(exc))
            sources[(lon, lat)] = source

    if missing:
        messages.addWarningMessage(
            "{} av {} tile(s) finns inte publicerade (vanligt över öppet vatten). "
            "Fyndpunkter där kan inte samplas.".format(len(missing), len(tiles))
        )
    if not sources:
        raise ValueError(
            "Ingen av fyndpunkternas tiles finns publicerad för {} i valt dataset. "
            "Kontrollera år och dataset-version.".format(year)
        )

    # Steg 2, huvudtråden: georeferering ur landmaskerna och rad/kolumn per fynd.
    work = []
    for source in sources.values():
        source.load_georeference()
    for obs in observations:
        source = sources.get(obs.tile)
        if source is None:
            obs.status = STATUS_NOTILE
            continue
        point = _project_geometry(obs.shape, source.sr).centroid
        obs.row, obs.col = source.rowcol(point.X, point.Y)
        work.append((obs, source, point.X, point.Y))

    if not work:
        raise ValueError(
            "Ingen fyndpunkt ligger i en publicerad tile för {} i valt dataset.".format(year)
        )

    # Steg 3, parallellt: två Range-anrop per fyndpunkt och radlöpa.
    arcpy.SetProgressor("step", "Hämtar embeddings för fyndpunkterna...", 0, len(work), 1)
    done = 0
    failures = []
    try:
        with concurrent.futures.ThreadPoolExecutor(max_workers=_WORKERS) as pool:
            futures = {
                pool.submit(_sample_one, source, obs, x, y, radius): obs
                for obs, source, x, y in work
            }
            for future in concurrent.futures.as_completed(futures):
                obs = futures[future]
                try:
                    vector = future.result()
                except _TileMissing:
                    vector = None
                except (urllib.error.URLError, OSError) as exc:
                    failures.append(exc)
                    vector = None
                if vector is not None:
                    obs.vector = vector
                    obs.status = STATUS_USED
                done += 1
                arcpy.SetProgressorPosition(done)
    finally:
        arcpy.ResetProgressor()

    if failures:
        raise ValueError(
            "{} fyndpunkt(er) kunde inte hämtas. {}".format(
                len(failures), _http_error_msg(failures[0]))
        )

    used = [obs for obs in observations if obs.status == STATUS_USED]
    if not used:
        raise ValueError(
            "Ingen fyndpunkt gav en giltig embedding-vektor. Ligger punkterna "
            "inom Tesseras täckning för valt år och dataset?"
        )

    no_data = sum(1 for obs in observations if obs.status == STATUS_NODATA)
    if no_data:
        messages.addWarningMessage(
            "{} fyndpunkt(er) låg på NoData i Tessera och hoppades över.".format(no_data)
        )
    if any(obs.clipped for obs in observations):
        messages.addWarningMessage(
            "Samplingsfönstret för minst en fyndpunkt nådde utanför sin tile och "
            "klipptes. Delen i grantile hämtades inte."
        )
    messages.addMessage(
        "{} av {} fyndpunkter fick en embedding-vektor.".format(
            len(used), len(observations))
    )
    return used


# =============================================================================
# Referensvektorer
# =============================================================================

def _unit_rows(matrix):
    """Radvis normerade vektorer. Nollvektorer lämnas orörda och fångas av
    anroparen."""
    norms = np.linalg.norm(matrix, axis=1)
    safe = np.where(norms == 0, 1.0, norms)
    return (matrix / safe[:, np.newaxis]).astype(np.float32), norms


def _neighbour_scores(units, k):
    """
    Den k:e största kosinuslikheten till någon annan fyndpunkt, per fyndpunkt.

    Det här är avsiktligt inte avstånd till tyngdpunkten. En art som använder
    två habitat är tvåtoppig i embedding-rummet, och en tyngdpunktsregel skulle
    kasta det mindre habitatet i sin helhet. Ett grannskapsmått frågar i stället
    om punkten har sällskap: en punkt inuti en verklig ansamling har nära
    grannar hur liten ansamlingen än är, medan ett felplacerat fynd ute till
    havs eller på ett hustak inte har några alls.

    Matrisen beräknas blockvis så att några tusen fyndpunkter ryms i minnet.
    """
    n = units.shape[0]
    scores = np.empty(n, dtype=np.float32)
    block = max(1, int(_KNN_BUDGET_BYTES / max(n * 4, 1)))
    for start in range(0, n, block):
        stop = min(n, start + block)
        sims = (units[start:stop] @ units.T).astype(np.float32)
        for offset in range(stop - start):
            sims[offset, start + offset] = -np.inf   # inte sig själv
        scores[start:stop] = np.partition(sims, -k, axis=1)[:, -k]
    return scores


def _drop_outliers(used, drop_pct, messages):
    """Uteslut de drop_pct procent fyndpunkter som har glesast grannskap."""
    units, norms = _unit_rows(np.vstack([obs.vector for obs in used]))
    if (norms == 0).any():
        raise ValueError("En fyndpunkts embedding-vektor är nollvektorn.")

    n = len(used)
    if n < 4:
        messages.addWarningMessage(
            "Färre än fyra fyndpunkter — utrensningen av avvikande fynd hoppades över."
        )
        return used

    k = min(max(3, int(math.ceil(0.05 * n))), 20, n - 1)
    scores = _neighbour_scores(units, k)
    for obs, score in zip(used, scores):
        obs.neighbour = float(score)

    n_drop = int(math.floor(n * drop_pct / 100.0))
    n_drop = min(n_drop, n - 1)
    if n_drop <= 0:
        messages.addMessage(
            "Grannskapsmått (k={}): inga fyndpunkter uteslöts.".format(k)
        )
        return used

    order = np.argsort(scores, kind="stable")
    dropped = set(int(i) for i in order[:n_drop])
    for position, obs in enumerate(used):
        if position in dropped:
            obs.status = STATUS_DROPPED
    kept = [obs for position, obs in enumerate(used) if position not in dropped]
    messages.addMessage(
        "Grannskapsmått (k={}): {} av {} fyndpunkter uteslöts som avvikande "
        "(likhet till {}:e grannen under {:.3f}).".format(
            k, n_drop, n, k, float(scores[order[n_drop]]))
    )
    return kept


def _reference_matrix(vectors, normalize, aggregation):
    """
    Referensmatrisen som varje pixel jämförs mot, en rad per referens.

    Med "största likhet" behålls samtliga fyndvektorer och varje pixel tar sin
    bästa träff. Med "medelvektor" kollapsar de till en rad: medelvärdet av
    likheten mot N normerade vektorer är exakt skalärprodukten mot medelvärdet
    av dem, så alternativet kostar en skalärprodukt i stället för N. Det är
    samma sammanvägning som "Tessera similarity search" gör med flera
    referensobjekt.
    """
    matrix = np.vstack(vectors).astype(np.float32)
    if normalize:
        matrix, norms = _unit_rows(matrix)
        if (norms == 0).any():
            raise ValueError("En referensvektor är nollvektorn och kan inte normaliseras.")
    if aggregation == AGG_MEAN:
        matrix = matrix.mean(axis=0, keepdims=True).astype(np.float32)
    return matrix


# =============================================================================
# Likhetsraster
# =============================================================================

def _similarity_raster(raster, reference, normalize, messages):
    """
    Likheten mot referensmatrisen för varje pixel: skalärprodukten mot varje rad,
    och den största av dem. Med en enda rad är det just den radens likhet.

    Indata läses radvis i block; blockets höjd sätts av både bandantalet och
    antalet referensvektorer, eftersom likhetsmatrisen (N x pixlar) kan bli
    större än det inlästa blocket när många fyndpunkter används.

    Anpassad från _similarity_raster i "Tessera similarity search".
    """
    extent = raster.extent
    width, height = raster.width, raster.height
    cell_w, cell_h = raster.meanCellWidth, raster.meanCellHeight
    n_bands = raster.bandCount
    n_refs = reference.shape[0]

    out = np.full((height, width), np.nan, dtype=np.float32)

    bytes_per_row = max((n_bands + n_refs) * width * 4, 1)
    block_rows = max(1, min(height, int(_BLOCK_BUDGET_BYTES / bytes_per_row)))

    arcpy.SetProgressor("step", "Beräknar likhet...", 0, height, block_rows)
    try:
        row = 0
        while row < height:
            rows = min(block_rows, height - row)
            y_bottom = extent.YMax - (row + rows) * cell_h
            origin = arcpy.Point(extent.XMin, y_bottom)

            block = arcpy.RasterToNumPyArray(
                raster, origin, width, rows, nodata_to_value=np.nan)
            if n_bands == 1:
                block = block[np.newaxis, :, :]
            block = block.astype(np.float32)

            invalid = np.isnan(block).any(axis=0)
            filled = np.where(np.isnan(block), np.float32(0.0), block)
            flat = filled.reshape(n_bands, rows * width)

            if normalize:
                norm = np.sqrt((flat ** 2).sum(axis=0))
                zero = norm == 0
                invalid = invalid | zero.reshape(rows, width)
                norm = np.where(zero, np.float32(1.0), norm)

            dots = reference @ flat
            if normalize:
                dots /= norm
            best = dots.max(axis=0).reshape(rows, width)
            best[invalid] = np.nan

            # Radblocket motsvarar rader [row, row + rows) uppifrån, samma
            # ordning som RasterToNumPyArray läser och NumPyArrayToRaster
            # förväntar sig vid skrivningen.
            out[row:row + rows, :] = best
            row += rows
            arcpy.SetProgressorPosition(row)
    finally:
        arcpy.ResetProgressor()

    return out


def _observation_scores(reference, vectors, normalize, aggregation):
    """
    Likhetsvärdet varje fyndpunkt själv skulle få, beräknat utan sitt eget
    bidrag till referensen.

    Utan den utelämningen blir svaret meningslöst: med "största likhet" matchar
    varje fyndpunkt sig själv och får 1,0. Värdena används för att kalibrera
    tröskeln — "90 procent av de kända fynden ligger över X".
    """
    matrix = np.vstack(vectors).astype(np.float32)
    if normalize:
        matrix, _ = _unit_rows(matrix)
    sims = (reference @ matrix.T).astype(np.float32)   # (referenser, fynd)

    n = matrix.shape[0]
    if aggregation == AGG_MEAN:
        # En enda referensrad, medelvärdet av alla fynd. Fyndets eget bidrag
        # räknas bort ur medelvärdet.
        own = (matrix * matrix).sum(axis=1) if not normalize else np.ones(n, dtype=np.float32)
        if n < 2:
            return sims[0]
        return ((sims[0] * n) - own) / (n - 1)

    if n < 2:
        return sims.max(axis=0)
    sims = sims.copy()
    sims[np.arange(n), np.arange(n)] = -np.inf
    return sims.max(axis=0)


def _save_raster(array, raster, raster_sr, out_path, messages):
    previous_sr = arcpy.env.outputCoordinateSystem
    arcpy.env.outputCoordinateSystem = raster_sr
    try:
        result = arcpy.NumPyArrayToRaster(
            array, arcpy.Point(raster.extent.XMin, raster.extent.YMin),
            raster.meanCellWidth, raster.meanCellHeight,
            value_to_nodata=np.nan,
        )
        result.save(out_path)
    finally:
        arcpy.env.outputCoordinateSystem = previous_sr

    try:
        arcpy.management.CalculateStatistics(out_path)
    except Exception as exc:
        messages.addWarningMessage("Kunde inte beräkna statistik för rastret: {}".format(exc))
    return out_path


def _threshold_polygons(sim_array, raster, raster_sr, threshold, out_path,
                        scratch_dir, messages):
    mask = np.where(np.isfinite(sim_array) & (sim_array >= threshold), 1, 0).astype(np.uint8)
    if not mask.any():
        messages.addWarningMessage(
            "Inga pixlar nådde tröskelvärdet {:.4f} — inga polygoner skapades.".format(threshold)
        )
        return None

    mask_path = os.path.join(scratch_dir, "habitat_mask.tif")
    previous_sr = arcpy.env.outputCoordinateSystem
    arcpy.env.outputCoordinateSystem = raster_sr
    try:
        result = arcpy.NumPyArrayToRaster(
            mask, arcpy.Point(raster.extent.XMin, raster.extent.YMin),
            raster.meanCellWidth, raster.meanCellHeight,
            value_to_nodata=0,
        )
        result.save(mask_path)
    finally:
        arcpy.env.outputCoordinateSystem = previous_sr

    arcpy.conversion.RasterToPolygon(mask_path, out_path, "SIMPLIFY", "VALUE")
    return out_path


# =============================================================================
# Utdatapunkter
# =============================================================================

def _write_points(observations, out_path, channels, keep_vectors, messages):
    """
    Fyndpunkterna med sitt likhetsvärde och sin status.

    Källans övriga attribut följer inte med — fyndets ObjectID sparas i
    kalla_oid, så en koppling tillbaka är en join bort. Det är avsiktligt:
    ett CopyFeatures och en matchning på radordning hade varit skörare än en
    uttrycklig nyckel.
    """
    workspace, name = os.path.dirname(out_path), os.path.basename(out_path)
    sr = observations[0].shape.spatialReference

    arcpy.management.CreateFeatureclass(
        workspace, name, "POINT", spatial_reference=sr)

    fields = [
        ["kalla_oid", "LONG", "Källans ObjectID"],
        ["grid", "TEXT", "Tessera-tile", 32],
        ["status", "TEXT", "Status", 24],
        ["likhet", "DOUBLE", "Likhet (utan eget bidrag)"],
        ["grannlikhet", "DOUBLE", "Likhet till k:e grannen"],
    ]
    channel_names = []
    if keep_vectors:
        for channel in channels:
            field = "kanal_{:03d}".format(channel + 1)
            channel_names.append(field)
            fields.append([field, "FLOAT", "Tessera-kanal {}".format(channel + 1)])
    arcpy.management.AddFields(out_path, fields)

    columns = ["SHAPE@", "kalla_oid", "grid", "status", "likhet", "grannlikhet"]
    columns.extend(channel_names)
    with arcpy.da.InsertCursor(out_path, columns) as cursor:
        for obs in observations:
            row = [
                obs.shape,
                obs.oid,
                _grid_name(*obs.tile),
                obs.status,
                obs.score,
                obs.neighbour,
            ]
            if keep_vectors:
                if obs.vector is None:
                    row.extend([None] * len(channel_names))
                else:
                    row.extend(float(v) for v in obs.vector[channels])
            cursor.insertRow(row)

    messages.addMessage("Skrev {}".format(out_path))
    return out_path


# =============================================================================
# Projekt- och standardvärden
# =============================================================================

def _default_gdb():
    try:
        aprx = arcpy.mp.ArcGISProject("CURRENT")
        if aprx.defaultGeodatabase:
            return aprx.defaultGeodatabase
    except Exception:
        pass
    workspace = arcpy.env.workspace
    if workspace and str(workspace).lower().endswith(".gdb"):
        return workspace
    return None


def _default_cache_dir():
    """
    Standardmapp för landmasker, skapad om den saknas.

    Inte tempfile.gettempdir(): inne i Pro pekar den på en egen mapp per session
    (ArcGISProTemp<nnnn>) som ofta blir kvar när Pro stängs. En sådan sökväg
    finns inte förrän någon skapar den, vilket får DEFolder-parametern att falla
    på ERROR 000732 redan när dialogen öppnas, och namnet byts vid varje omstart
    så att ingenting återanvänds.
    """
    base = os.environ.get("LOCALAPPDATA") or tempfile.gettempdir()
    path = os.path.join(base, _CACHE_DIRNAME)
    try:
        os.makedirs(path, exist_ok=True)
    except OSError:
        pass
    return path


def _looks_synced(path):
    lowered = (path or "").lower()
    return any(hint in lowered for hint in _SYNC_HINTS)


def _year_in_name(text):
    """Årtalet i ett rasternamn, t.ex. 'tessera_2024' -> '2024'."""
    for match in re.findall(r"(20\d{2})", text or ""):
        if match in YEARS:
            return match
    return None


def _add_to_map(outputs, messages):
    try:
        aprx = arcpy.mp.ArcGISProject("CURRENT")
        map_obj = aprx.activeMap
        if map_obj is None:
            maps = aprx.listMaps()
            map_obj = maps[0] if maps else None
    except Exception:
        map_obj = None
    if map_obj is None:
        messages.addWarningMessage("Ingen aktiv karta — resultatet lades inte till.")
        return
    for path in outputs:
        try:
            map_obj.addDataFromPath(path)
        except Exception as exc:
            messages.addWarningMessage("  Kunde inte lägga till {}: {}".format(path, exc))


# =============================================================================
# Kanalangivelse
# =============================================================================

def _parse_channels(text, count=N_CHANNELS):
    """
    Tolka en kanalangivelse som "1-16,64" till nollbaserade index. Tom sträng
    ger alla kanaler. Kanalerna numreras 1-128, som banden i dialogen.
    """
    text = (text or "").strip()
    if not text:
        return list(range(count))

    if not re.fullmatch(r"[0-9,\-\s]+", text):
        raise ValueError(
            "Ogiltig kanalangivelse: '{}'. Ange kanaler som t.ex. 1-16,64.".format(text)
        )

    indices = []
    for part in text.split(","):
        part = part.strip()
        if not part:
            continue
        match = re.fullmatch(r"(\d+)\s*-\s*(\d+)", part)
        if match:
            first, last = int(match.group(1)), int(match.group(2))
            if first > last:
                raise ValueError("Ogiltigt kanalintervall: '{}'.".format(part))
            values = range(first, last + 1)
        else:
            values = [int(part)]
        for value in values:
            if not 1 <= value <= count:
                raise ValueError(
                    "Kanal {} finns inte — Tessera har kanal 1-{}.".format(value, count)
                )
            if value - 1 not in indices:
                indices.append(value - 1)

    if not indices:
        raise ValueError("Ingen giltig kanalangivelse.")
    return sorted(indices)


def _channels_for_raster(text, band_count):
    """
    Vilka Tessera-kanaler sökrastrets band motsvarar.

    Bandkopplingen går inte att läsa ur rastret. "Tessera embeddings to GDB"
    låter användaren spara ett urval av kanalerna, så band 3 i ett 16-bands
    raster kan vara kanal 64. De samplade punktvektorerna är alltid 128 långa
    och måste skäras ned till samma kanaler, i samma ordning.
    """
    channels = _parse_channels(text)
    if len(channels) != band_count:
        raise ValueError(
            "Sökrastret har {} band men kanalangivelsen räknar upp {} kanaler. "
            "Ange vilka Tessera-kanaler rastrets band motsvarar, t.ex. 1-{}.".format(
                band_count, len(channels), band_count)
        )
    return channels


# =============================================================================
# Toolbox
# =============================================================================

class Toolbox:
    def __init__(self):
        self.label = "Tessera habitat"
        self.alias = "tesserahabitat"
        self.tools = [HabitatPrediktion]


class HabitatPrediktion:
    def __init__(self):
        self.label = "Habitatprediktion från fyndpunkter"
        self.description = (
            "Hämtar Tessera-embeddingen i varje fyndpunkt för en art och räknar ut "
            "hur lik varje pixel i ett sökraster är dessa kända lokaler. Resultatet "
            "är ett enbandsraster där höga värden betyder att platsen ser ut som de "
            "platser arten redan hittats på.\n\n"
            "Fyndpunkterna hämtas med byte-intervall direkt ur Tesseras publicerade "
            "filer, så ingen tile behöver laddas ned — kostnaden är ett par små "
            "anrop per fyndpunkt oavsett hur utspridda de är.\n\n"
            "Standardvalet jämför varje pixel med samtliga fyndpunkter och behåller "
            "den bästa träffen, vilket bevarar flera habitat. Fyndens egna "
            "likhetsvärden redovisas i loggen som stöd för att välja tröskel."
        )
        self.canRunInBackground = False
        self._raster_memo = ""
        self._poly_memo = ""
        self._year_memo = DEFAULT_YEAR

    # ── Parametrar ────────────────────────────────────────────────────────────

    def getParameterInfo(self):
        p_points = arcpy.Parameter(
            displayName="Fyndpunkter (observationer av en art)",
            name="obs_points", datatype="GPFeatureLayer",
            parameterType="Required", direction="Input",
        )
        p_points.filter.list = ["Point"]

        p_raster = arcpy.Parameter(
            displayName="Embedding-raster att söka i (sökområdet)",
            name="search_raster", datatype="DERasterDataset",
            parameterType="Required", direction="Input",
        )

        p_channels = arcpy.Parameter(
            displayName="Tessera-kanaler i sökrastret, t.ex. 1-16,64 (tomt = 1-128)",
            name="raster_channels", datatype="GPString",
            parameterType="Optional", direction="Input",
        )

        p_year = arcpy.Parameter(
            displayName="År",
            name="year", datatype="GPString",
            parameterType="Required", direction="Input", category=CAT_TESSERA,
        )
        p_year.filter.type = "ValueList"
        p_year.filter.list = YEARS
        p_year.value = DEFAULT_YEAR

        p_dataset = arcpy.Parameter(
            displayName="Dataset-version",
            name="dataset", datatype="GPString",
            parameterType="Required", direction="Input", category=CAT_TESSERA,
        )
        p_dataset.filter.type = "ValueList"
        p_dataset.filter.list = list(DATASETS)
        p_dataset.value = DEFAULT_DATASET

        p_radius = arcpy.Parameter(
            displayName="Samplingsradie runt varje fyndpunkt (m, 0 = en pixel)",
            name="sample_radius", datatype="GPDouble",
            parameterType="Optional", direction="Input", category=CAT_TESSERA,
        )
        p_radius.value = 0.0

        p_cache = arcpy.Parameter(
            displayName="Cache-mapp för landmasker",
            name="cache_dir", datatype="DEFolder",
            parameterType="Optional", direction="Input", category=CAT_TESSERA,
        )
        p_cache.value = _default_cache_dir()

        p_agg = arcpy.Parameter(
            displayName="Sammanvägning av fyndpunkterna",
            name="aggregation", datatype="GPString",
            parameterType="Required", direction="Input", category=CAT_REF,
        )
        p_agg.filter.type = "ValueList"
        p_agg.filter.list = [AGG_MAX, AGG_MEAN]
        p_agg.value = AGG_MAX

        p_drop = arcpy.Parameter(
            displayName="Uteslut avvikande fyndpunkter (% med glesast grannskap)",
            name="drop_pct", datatype="GPDouble",
            parameterType="Optional", direction="Input", category=CAT_REF,
        )
        p_drop.value = 0.0

        p_normalize = arcpy.Parameter(
            displayName="Normalisera vektorer (kosinuslikhet, -1 till 1)",
            name="normalize", datatype="GPBoolean",
            parameterType="Optional", direction="Input", category=CAT_REF,
        )
        p_normalize.value = True

        p_out_raster = arcpy.Parameter(
            displayName="Utdata: likhetsraster",
            name="out_raster", datatype="DERasterDataset",
            parameterType="Required", direction="Output",
        )
        gdb = _default_gdb()
        if gdb:
            # Startvärdet skrivs också till memot, annars ser updateParameters
            # det som ett namn användaren själv skrivit och låter bli att
            # föreslå <raster>_habitat.
            self._raster_memo = os.path.join(gdb, "habitat")
            p_out_raster.value = self._raster_memo

        p_out_points = arcpy.Parameter(
            displayName="Utdata: fyndpunkter med likhetsvärde",
            name="out_points", datatype="DEFeatureClass",
            parameterType="Optional", direction="Output", category=CAT_OUT,
        )

        p_keep = arcpy.Parameter(
            displayName="Spara embedding-värdena som fält i utdatapunkterna",
            name="keep_vectors", datatype="GPBoolean",
            parameterType="Optional", direction="Input", category=CAT_OUT,
        )
        p_keep.value = False

        p_pct = arcpy.Parameter(
            displayName="Tröskel som percentil av fyndens likhet, t.ex. 10 (tomt = hoppa över)",
            name="pct_threshold", datatype="GPDouble",
            parameterType="Optional", direction="Input", category=CAT_OUT,
        )

        p_threshold = arcpy.Parameter(
            displayName="Eller absolut tröskelvärde",
            name="threshold", datatype="GPDouble",
            parameterType="Optional", direction="Input", category=CAT_OUT,
        )

        p_out_polygons = arcpy.Parameter(
            displayName="Utdata: polygoner över tröskelvärdet",
            name="out_polygons", datatype="DEFeatureClass",
            parameterType="Optional", direction="Output", category=CAT_OUT,
        )

        p_overwrite = arcpy.Parameter(
            displayName="Skriv över befintlig utdata",
            name="overwrite", datatype="GPBoolean",
            parameterType="Optional", direction="Input", category=CAT_OUT,
        )
        p_overwrite.value = True

        p_add = arcpy.Parameter(
            displayName="Lägg till resultatet i kartan",
            name="add_to_map", datatype="GPBoolean",
            parameterType="Optional", direction="Input", category=CAT_MAP,
        )
        p_add.value = True

        return [p_points, p_raster, p_channels, p_year, p_dataset, p_radius,
                p_cache, p_agg, p_drop, p_normalize, p_out_raster, p_out_points,
                p_keep, p_pct, p_threshold, p_out_polygons, p_overwrite, p_add]

    def isLicensed(self):
        return True

    # ── Dialog ────────────────────────────────────────────────────────────────

    def updateParameters(self, parameters):
        (_p_points, p_raster, _p_channels, p_year, _p_dataset, _p_radius,
         _p_cache, _p_agg, _p_drop, _p_normalize, p_out_raster, p_out_points,
         p_keep, p_pct, p_threshold, p_out_polygons, _p_overwrite,
         _p_add) = parameters

        if p_raster.value is not None:
            in_name = os.path.basename(str(p_raster.valueAsText)).rsplit(".", 1)[0]

            # Namnförslaget följer indatarastret tills användaren skrivit ett
            # eget. Jämförelsen görs mot det senast föreslagna namnet, inte mot
            # altered-flaggan, som också sätts när koden själv skriver värdet.
            gdb = _default_gdb()
            suggestion = os.path.join(gdb, "{}_habitat".format(in_name)) if gdb else ""
            if suggestion and (p_out_raster.valueAsText or "").strip() in ("", self._raster_memo):
                p_out_raster.value = suggestion
            self._raster_memo = suggestion

            # Året gissas ur rasternamnet: "Tessera embeddings to GDB" döper
            # sina raster till tessera_<år>, och punktsamplingen måste hämta
            # samma år som rastret innehåller.
            found = _year_in_name(in_name)
            if found and (p_year.valueAsText or "") == self._year_memo:
                p_year.value = found
            self._year_memo = p_year.valueAsText or DEFAULT_YEAR

        p_keep.enabled = bool((p_out_points.valueAsText or "").strip())

        has_threshold = p_pct.value is not None or p_threshold.value is not None
        p_out_polygons.enabled = has_threshold
        if has_threshold:
            out_text = (p_out_raster.valueAsText or "").strip()
            if out_text:
                suggestion = out_text + "_omraden"
                if (p_out_polygons.valueAsText or "").strip() in ("", self._poly_memo):
                    p_out_polygons.value = suggestion
                self._poly_memo = suggestion

    def updateMessages(self, parameters):
        (p_points, p_raster, p_channels, p_year, p_dataset, p_radius,
         p_cache, _p_agg, p_drop, p_normalize, _p_out_raster, _p_out_points,
         _p_keep, p_pct, p_threshold, p_out_polygons, _p_overwrite,
         _p_add) = parameters

        if p_points.value is not None:
            try:
                if arcpy.Describe(p_points.value).shapeType != "Point":
                    p_points.setErrorMessage("Fyndlagret måste innehålla punkter.")
            except Exception:
                pass

        band_count = None
        if p_raster.value is not None:
            try:
                band_count = int(arcpy.Describe(p_raster.value).bandCount)
            except Exception:
                band_count = None

        if band_count is not None:
            try:
                _channels_for_raster(p_channels.valueAsText, band_count)
            except ValueError as exc:
                p_channels.setErrorMessage(str(exc))

            in_name = os.path.basename(str(p_raster.valueAsText)).rsplit(".", 1)[0]
            found = _year_in_name(in_name)
            if found and found != (p_year.valueAsText or ""):
                p_year.setWarningMessage(
                    "Rastret heter {} men punkterna samplas för {}. En blandning av "
                    "årgångar ger fel svar utan att fela.".format(in_name, p_year.valueAsText)
                )

        if p_dataset.valueAsText and p_dataset.valueAsText != DEFAULT_DATASET:
            p_dataset.setWarningMessage(
                "Olika Tessera-versioner har olika kanalrum. Sökrastret måste vara "
                "hämtat ur samma version som punkterna samplas ur."
            )

        if p_radius.value is not None:
            if p_radius.value < 0:
                p_radius.setErrorMessage("Samplingsradien kan inte vara negativ.")
            elif p_radius.value > _MAX_RADIUS_M:
                p_radius.setErrorMessage(
                    "Samplingsradien får vara högst {:.0f} m.".format(_MAX_RADIUS_M)
                )
            elif p_radius.value > 100:
                p_radius.setWarningMessage(
                    "En stor radie blandar in omgivande habitat i fyndets vektor."
                )

        if p_drop.value is not None and not (0 <= p_drop.value < 100):
            p_drop.setErrorMessage("Andelen måste ligga mellan 0 och 100.")

        if p_pct.value is not None and not (0 <= p_pct.value <= 100):
            p_pct.setErrorMessage("Percentilen måste ligga mellan 0 och 100.")

        if p_pct.value is not None and p_threshold.value is not None:
            message = "Ange antingen en percentil eller ett absolut tröskelvärde, inte båda."
            p_pct.setErrorMessage(message)
            p_threshold.setErrorMessage(message)

        if (p_threshold.value is not None
                and bool(p_normalize.value if p_normalize.value is not None else True)
                and not (-1.0 <= p_threshold.value <= 1.0)):
            p_threshold.setWarningMessage(
                "Med normaliserade vektorer ligger likheten mellan -1 och 1."
            )

        if ((p_pct.value is not None or p_threshold.value is not None)
                and not (p_out_polygons.valueAsText or "").strip()):
            p_out_polygons.setErrorMessage(
                "Ange var polygonerna över tröskelvärdet ska sparas."
            )

        if p_cache.valueAsText and _looks_synced(p_cache.valueAsText):
            p_cache.setWarningMessage(
                "Mappen ser ut att synkas till molnet. Välj hellre en lokal mapp."
            )

    # ── Körning ───────────────────────────────────────────────────────────────

    def execute(self, parameters, messages):
        try:
            _run(
                points=parameters[0].valueAsText,
                raster_path=parameters[1].valueAsText,
                channels_text=parameters[2].valueAsText,
                year=parameters[3].valueAsText,
                dataset=parameters[4].valueAsText,
                radius=parameters[5].value or 0.0,
                cache_dir=parameters[6].valueAsText,
                aggregation=parameters[7].valueAsText,
                drop_pct=parameters[8].value or 0.0,
                normalize=bool(parameters[9].value) if parameters[9].value is not None else True,
                out_raster_path=parameters[10].valueAsText,
                out_points_path=parameters[11].valueAsText,
                keep_vectors=bool(parameters[12].value),
                pct_threshold=parameters[13].value,
                threshold=parameters[14].value,
                out_polygons_path=parameters[15].valueAsText,
                overwrite=bool(parameters[16].value) if parameters[16].value is not None else True,
                add_to_map=bool(parameters[17].value) if parameters[17].value is not None else True,
                messages=messages,
            )
        except ValueError as exc:
            messages.addErrorMessage(str(exc))
            raise arcpy.ExecuteError

    def postExecute(self, parameters):
        return


# =============================================================================
# Körningens innehåll (separat funktion — går att testa utanför Pro)
# =============================================================================

def _check_output(path, overwrite):
    if arcpy.Exists(path):
        if not overwrite:
            raise ValueError(
                "{} finns redan. Kryssa i 'Skriv över befintlig utdata' eller välj "
                "ett annat namn.".format(path)
            )
        arcpy.management.Delete(path)


def _run(points, raster_path, channels_text, year, dataset, radius, cache_dir,
         aggregation, drop_pct, normalize, out_raster_path, out_points_path,
         keep_vectors, pct_threshold, threshold, out_polygons_path, overwrite,
         add_to_map, messages):
    """Hela körningen. Returnerar listan med skapade dataset."""

    if not arcpy.Exists(points):
        raise ValueError("Fyndlagret {} finns inte.".format(points))
    if not arcpy.Exists(raster_path):
        raise ValueError("Sökrastret {} finns inte.".format(raster_path))
    if not out_raster_path:
        raise ValueError("Ange var likhetsrastret ska sparas.")

    if dataset not in DATASETS:
        raise ValueError("Okänd dataset-version: {}.".format(dataset))

    # Cache-mappen måste vara bestämd innan klienten skapas: den avgör var
    # geotessera lägger landmaskerna.
    cache_dir = cache_dir or _default_cache_dir()
    if _looks_synced(cache_dir):
        messages.addWarningMessage(
            "Cache-mappen {} ser ut att synkas till molnet.".format(cache_dir)
        )
    os.makedirs(cache_dir, exist_ok=True)
    client = _client(dataset, cache_dir, messages)

    if year not in YEARS:
        raise ValueError("Okänt år: {}.".format(year))
    aggregation = aggregation or AGG_MAX
    if aggregation not in (AGG_MAX, AGG_MEAN):
        raise ValueError("Okänd sammanvägning: {}.".format(aggregation))

    if pct_threshold is not None and threshold is not None:
        raise ValueError(
            "Ange antingen en percentil eller ett absolut tröskelvärde, inte båda."
        )
    if (pct_threshold is not None or threshold is not None) and not out_polygons_path:
        raise ValueError("Ange var polygonerna över tröskelvärdet ska sparas.")

    _check_output(out_raster_path, overwrite)
    if out_points_path:
        _check_output(out_points_path, overwrite)
    if out_polygons_path and (pct_threshold is not None or threshold is not None):
        _check_output(out_polygons_path, overwrite)

    raster = arcpy.Raster(raster_path)
    raster_sr = raster.spatialReference
    if not _sr_is_valid(raster_sr):
        raise ValueError("Sökrastret saknar koordinatsystem.")

    channels = _channels_for_raster(channels_text, raster.bandCount)
    messages.addMessage(
        "Sökrastrets {} band tolkas som Tessera-kanal {}.".format(
            raster.bandCount,
            "1-{}".format(len(channels)) if channels == list(range(len(channels)))
            else ", ".join(str(c + 1) for c in channels[:8]) + ("..." if len(channels) > 8 else ""))
    )
    messages.addMessage(
        "Fyndpunkterna samplas ur Tessera {}, år {}.".format(dataset, year)
    )

    if raster.width * raster.height > _MAX_SANE_PIXELS:
        messages.addWarningMessage(
            "Sökrastret har {} x {} pixlar. Likhetsrastret hålls helt i minnet — "
            "klipp indatarastret till ett mindre område om du får minnesfel.".format(
                raster.width, raster.height)
        )

    observations = _read_observations(points, messages)
    used = _sample_observations(
        observations, client, year, radius or 0.0, messages)

    if drop_pct and drop_pct > 0:
        used = _drop_outliers(used, drop_pct, messages)

    # Fyndvektorerna skärs ned till samma kanaler som sökrastret innehåller.
    vectors = [obs.vector[channels] for obs in used]
    reference = _reference_matrix(vectors, normalize, aggregation)
    messages.addMessage(
        "Referens: {} ({} vektor(er), {} kanaler).".format(
            "största likhet" if aggregation == AGG_MAX else "medelvektor",
            reference.shape[0], reference.shape[1])
    )

    scratch_dir = tempfile.mkdtemp(prefix=_SCRATCH_DIRNAME + "_")
    outputs = []
    try:
        messages.addMessage("Beräknar likhet mot hela sökrastret...")
        sim_array = _similarity_raster(raster, reference, normalize, messages)

        valid = np.isfinite(sim_array)
        if not valid.any():
            raise ValueError("Alla pixlar blev NoData — inget resultat att skriva.")
        messages.addMessage(
            "Likhet ({}): min {:.3f}, max {:.3f}.".format(
                "kosinus" if normalize else "oskalad skalärprodukt",
                float(sim_array[valid].min()), float(sim_array[valid].max()))
        )

        # Fyndens egna värden, utan eget bidrag till referensen. De är
        # kalibreringen: "90 procent av de kända fynden ligger över X".
        scores = _observation_scores(reference, vectors, normalize, aggregation)
        for obs, score in zip(used, scores):
            obs.score = float(score)
        percentiles = np.percentile(scores, [0, 10, 25, 50])
        messages.addMessage(
            "Fyndpunkternas egen likhet (utan eget bidrag): min {:.3f}, "
            "10:e percentilen {:.3f}, 25:e {:.3f}, median {:.3f}.".format(*percentiles)
        )

        cut = threshold
        if pct_threshold is not None:
            cut = float(np.percentile(scores, pct_threshold))
            messages.addMessage(
                "Tröskeln sätts till {:.4f} — {:.0f} % av fyndpunkterna ligger "
                "under det värdet.".format(cut, pct_threshold)
            )

        _save_raster(sim_array, raster, raster_sr, out_raster_path, messages)
        messages.addMessage("Skrev {}".format(out_raster_path))
        outputs.append(out_raster_path)

        if cut is not None:
            polygons = _threshold_polygons(
                sim_array, raster, raster_sr, cut, out_polygons_path,
                scratch_dir, messages)
            if polygons:
                messages.addMessage("Skrev {}".format(polygons))
                outputs.append(polygons)
    finally:
        shutil.rmtree(scratch_dir, ignore_errors=True)

    if out_points_path:
        _write_points(observations, out_points_path, channels, keep_vectors, messages)
        outputs.append(out_points_path)

    if add_to_map:
        _add_to_map(outputs, messages)

    messages.addMessage("Klar!")
    return outputs
