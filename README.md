# UVNasalViz
Visualize and analyze data in simulated nasal cavities

Modified from UVGUI produced by ...
  (Reference here)
  
See also
- Elad et al. 2005
- others?

CONTENTS: 
| File name | Description |  
| --------- | ----------- |
| **caseDirs** | Sample OpenFoam case directories
| caseDirs/Elad_Test_Case | Approximates "geometric" domain from Elad et al. 2005. |
| **code** |  Matlab routines |
| code/calcAvgOverCycles.m | Calculate averages over several breath cycles | 
| code/calcPrefPhase.m | Calculate preferred phase  | 
| code/test_facestats.m | |
| code/test_findb.m | |
| code/test_plot_UV.m | |
| code/test_read_polyMesh.m | |
| code/test_WSS_overTime.m | Read WSS files for a sequence of time steps|
| **code/util_files** | Routines for reading files |
| code/util_files/facestats.m | |
| code/util_files/read_faces.m | |
| code/util_files/read_polyMesh.m | |
| code/util_files/read_uv.m | |
| code/util_files/read_wss_mag_OF_face.m | |
| **code/util_geom** | Routines for interacting with geometric entities; finding boundaries, areas, etc..
| code/util_geom/compAllWallBdy.m | Compute each piece of wall boundary separately, return in a cell structure|
| code/util_geom/faceArea.m | | 
| code/util_geom/findb.m | | 
| code/util_geom/mapQuantPointToFace.m | Map a quantity from points to faces | 
| code/util_geom/mapQuantFaceToPoint.m | Map a quantity from faces to points | 
| code/util_geom/maxSide.m | | 
