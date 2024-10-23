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
| caseDirs/Inthavong_Test_Case | Human nasal cavity from Inthavong et al. 2014, 2017 |
| caseDirs/ReadTecplot_GetRegions | To identify olfactory region of Inthavong mesh |
| **code** |  Matlab routines |
| code/calcAvgOverCycles.m | Calculate averages over several breath cycles | 
| code/calcPrefPhase.m | Calculate preferred phase  | 
| code/makePlotsCRCNSPoster.m | Plots for CRCNS poster, October 2022 |
| code/makePlotsCRCTalk.m | Plots for CRC talk, Aoril 2023 |
| code/test_facestats.m | |
| code/test_findb.m | |
| code/test_plot_UV.m | |
| code/test_read_polyMesh.m | |
| code/test_useOfBdyFile.m | Read in boundary file (which assigns each face to a physical region) |
| code/test_WSS_overTime.m | Read WSS files for a sequence of time steps|
| **code/util_files** | Routines for reading files |
| code/util_files/facestats.m | |
| code/util_files/read_faces.m | |
| code/util_files/read_polyMesh.m | |
| code/util_files/read_uv.m | |
| code/util_files/read_wss_mag_OF_face.m | |
| **code/util_geom** | Routines for interacting with geometric entities; finding boundaries, areas, etc..
| code/util_geom/compAllWallBdy.m | Compute each piece of wall boundary separately, return in a cell structure |
| code/util_geom/compCentroidsUV.m | Compute centroids of a list of UV faces  | 
| code/util_geom/cutByPlane.m | Cut a mesh by a plane; keep only one side | 
| code/util_geom/faceArea.m | | 
| code/util_geom/findb.m | | 
| code/util_geom/getplane.m | Find all intersections of boundary with a plane | 
| code/util_geom/mapQuantPointToFace.m | Map a quantity from points to faces | 
| code/util_geom/mapQuantFaceToPoint.m | Map a quantity from faces to points | 
| code/util_geom/maxSide.m | | 
| **code/util_plot** | Routines for plotting
| code/util_plot/makeRedBluePhaseCmap.m | Make a colormap for preferred phase (actually uses blue, violet, red, orange in that order) |
| code/util_plot/plotAllWallBdy.m | Plot the boundary of each wall piece |
