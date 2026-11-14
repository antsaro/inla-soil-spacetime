# Input data

The data are not included in this repository. Put the input file here as:

    data/parcels_with_soil_complete.gpkg

## Format

A single GeoPackage layer, **one row per field and season** (a field can have more than one row in a season when several sub-plots were sampled). Field polygons must be in a projected CRS in metres, because the neighbourhood graph uses boundary-to-boundary distances (the original data use EPSG:32628, UTM 28N).

Columns used by the code:

| Column | Type | Meaning |
|---|---|---|
| `ID.Field` | integer | field index, **1 to n with no gaps** (used as the graph node index) |
| `Year` | integer | season |
| `Bloc` | text | management block (fixed effect) |
| `irrg_ty` | text | irrigation zone, e.g. `GRV` (gravity) / `GAG` (drip) (fixed effect) |
| `geometry` | polygon | field boundary |

Soil columns (empty when not measured). Any of these that are present are modelled:

`pH`, `CN_Ratio`, `Electrical_Conductivity_mmho_cm`, `Exchangeable_Calcium_meq_100g`, `Exchangeable_Magnesium_meq_100g`, `Exchangeable_Sodium_meq_100g`, `Exchangeable_Potassium_meq_100g`, `Available_P2O5_mg_kg`, `Organic_Carbon_per_mille`, `Total_Nitrogen_per_mille`, `Organic_Matter_percent`, `Sum_of_Exchangeable_Bases_meq_100g`, `CaO_g_100g`, `Mg_Ca_Ratio`, `Mg_K_Ratio`, `Compaction_Index`

Other numeric columns that are not identifiers or coordinates are also picked up. Variables with fewer than `MIN_OBS` measurements (200 by default) are skipped.

Texture (used only through ILR coordinates):

`Clay_percent`, `Fine_Silt_percent`, `Coarse_Silt_percent`, `Fine_Sand_percent`, `Coarse_Sand_percent`, `Total_Silt_percent`, `Total_Sand_percent`

Values below 0.1 % are set to 0.1 % and each composition is closed to 100 % before the log-ratios are taken.
