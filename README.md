# Embedding habitat prediction

ArcGIS Pro Python toolbox that predicts where a species could occur, from points where it is
already known to occur. It looks up the Tessera embedding at each observation, then scores every
pixel of a search raster by how much it resembles those known sites. High values mean the place
looks like the places the species has been found.

The tile registry and the landmasks come from the `geotessera` package. The observations
themselves do not: they are read straight out of Tessera's published files with HTTP range
requests, one pixel at a time. A tile is about 90 MB and geotessera always fetches whole tiles,
while occurrence data for a species is typically spread across a whole region, so a tile per
observation is not workable. The files are uncompressed `.npy`, so a single pixel's 128 channel
vector is 128 contiguous bytes at a computable offset. Measured cost per observation is about
0.6 kB, plus one 13.6 kB landmask per tile that is cached and reused.

The similarity itself is the same dot product method as "Tessera similarity search" in this
project, which follows Google Earth Engine's similarity search tutorial.

https://developers.google.com/earth-engine/tutorials/community/satellite-embedding-05-similarity-search

## Requirements

- ArcGIS Pro 3.x. Developed and tested on 3.6 with Python 3.13. No Spatial Analyst or 3D
  Analyst extension is needed.
- The `geotessera` package, for the tile registry and the landmasks.
- Internet access to `s3.us-west-2.amazonaws.com`.

`geotessera` is not part of the default `arcgispro-py3` environment, and ArcGIS Pro does not
allow installing into it. Clone the environment first.

## Install

1. In ArcGIS Pro: Settings, Package Manager, clone the active environment. Name the clone
   something like `arcgispro-py3-personal` and make it active.
2. Install geotessera into the clone, then pin pyarrow back to the version Pro supports.
   geotessera pulls a newer pyarrow that fails to load once `arcpy` is imported, which breaks
   the registry:

   ```
   python -m pip install geotessera
   python -m pip install "pyarrow==20.0.0"
   ```

3. Clone or download this repo.
4. In ArcGIS Pro: Catalog, Toolboxes, Add Toolbox, select `TesseraHabitat.pyt`.
5. Open Tessera habitat, Habitatprediktion från fyndpunkter.

If the environment is wrong the tool stops with a message naming the active environment rather
than failing part way through.

## The tool dialog

The UI is in Swedish, matching a Swedish ArcGIS Pro install.

| Parameter | Default | Notes |
|---|---|---|
| Fyndpunkter | - | Point layer, observations of one species |
| Embedding-raster att söka i | - | The search area. A Tessera raster, e.g. from "Tessera embeddings to GDB" |
| Tessera-kanaler i sökrastret | all 128 | Which Tessera channels the raster's bands are. Required when the raster does not have 128 bands |
| År | guessed from the raster name | Must match the year the search raster holds |
| Dataset-version | v1 | Must match the version the search raster came from |
| Samplingsradie | 0 | Metres around each observation. 0 samples one pixel |
| Cache-mapp för landmasker | system temp | Reused across runs. Avoid cloud-synced folders |
| Sammanvägning av fyndpunkterna | största likhet | Or a single mean vector |
| Uteslut avvikande fyndpunkter | 0 % | Drops the observations with the sparsest neighbourhood |
| Normalisera vektorer | on | Cosine similarity (-1 to 1) instead of a raw dot product |
| Utdata: likhetsraster | `<raster>_habitat` in the project default gdb | |
| Utdata: fyndpunkter med likhetsvärde | - | The observations with their own score and status |
| Spara embedding-värdena som fält | off | Adds one field per channel to the point output |
| Tröskel som percentil av fyndens likhet | - | E.g. 10, meaning 90 percent of known occurrences score above the threshold |
| Eller absolut tröskelvärde | - | Use one or the other, not both |
| Utdata: polygoner över tröskelvärdet | `<likhetsraster>_omraden` | Required if a threshold is set |
| Skriv över befintlig utdata | on | |
| Lägg till resultatet i kartan | on | |

## Output

A single band float raster on the same grid as the search raster, holding the similarity at every
pixel. Pixels that are NoData in any band of the input are NoData in the output.

With a threshold set, a second output holds the polygons where the similarity is at or above that
value. The point output carries each observation's `status` (använd, utesluten, ingen tile, saknar
data), its `grannlikhet` (the neighbourhood statistic used for outlier removal) and its `likhet`.
Source attributes do not follow along; the source ObjectID is in `kalla_oid` so a join back is one
step.

## Combining several observations

A species usually uses more than one habitat, so occurrence data is multimodal in embedding space.
The default, `Största likhet mot någon fyndpunkt`, compares each pixel against every observation
and keeps the best match. Two distinct habitats both stay visible. The alternative,
`Medelvektor`, collapses the observations to one vector; it is faster and steadier when there are
very many observations, but a species that uses both open sand and old forest ends up represented
by a point midway between them that matches neither.

## Outlier removal

Optional, off by default. The rule is deliberately not distance to the centroid. With a bimodal
distribution the centroid sits between the two habitats, and the observations furthest from it are
exactly the smaller habitat, so a centroid rule deletes a real habitat wholesale.

Instead each observation is scored by its similarity to its k-th nearest other observation, and the
lowest percent are dropped. An observation inside a genuine cluster keeps close neighbours however
small its cluster is; a mislocated record out at sea or on a roof has none. k is chosen
automatically as 5 percent of the observation count, at least 3 and at most 20, and is printed in
the run log.

In a test with 32 observations in a dominant habitat, 4 in a minority habitat and 4 scattered junk
records, this rule dropped all 4 junk records and kept the minority habitat intact. A centroid rule
on the same data dropped the entire minority habitat and kept every junk record.

## Choosing a threshold

The run log reports what the observations themselves score, computed without their own contribution
to the reference. Without that leave-one-out step the number is meaningless: under
`Största likhet` every observation matches itself and scores 1.0.

```
Fyndpunkternas egen likhet (utan eget bidrag): min 0.620, 10:e percentilen 0.640, 25:e 0.805, median 0.918
```

Read as: set the threshold at 0.640 and 90 percent of the known occurrences would fall inside the
predicted area. Entering 10 under `Tröskel som percentil av fyndens likhet` does exactly that.

## Notes on the data

- The search raster and the sampled observations must come from the same Tessera year and the same
  dataset version. Different versions have different 128 channel feature spaces and a mixture gives
  a wrong answer without failing. The tool warns when the raster name names a different year.
- Band alignment cannot be read out of a raster. "Tessera embeddings to GDB" can save a subset of
  channels, so band 3 of a 16 band raster may be Tessera channel 64. Tell the tool which channels
  the raster holds; it defaults to 1-128 and refuses to run when the count does not match the band
  count.
- Tessera tiles are the UTM bounding box of a 0.1 degree cell, so they overspill that cell and
  neighbouring tiles overlap slightly. An observation is always sampled from the tile whose 0.1
  degree cell contains it. Near a tile edge that grid can be offset by a fraction of a cell from
  the search raster's grid, which is a sub-pixel difference, not an error.
- Tiles over open water are not published and return HTTP 404. Those observations are reported as
  `ingen tile` and skipped.
- A sampling radius larger than zero averages the cells whose centres fall within it. A window that
  reaches past the tile edge is clipped to the tile and the run log says so.
- The embedding is stored quantised as int8 with one scale factor per pixel, shared by all 128
  channels. The tool multiplies it out. With normalization on, that scale cancels.

## Source

- Project: https://geotessera.org/
- Library: https://github.com/ucam-eo/geotessera
- Project: https://geotessera.org/
- Method: https://developers.google.com/earth-engine/tutorials/community/satellite-embedding-05-similarity-search
