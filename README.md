# 360-degree-video-streaming
* This source code is an implementation in the MATLAB language of four adaptation algorithms for 360-degree video streaming. 
* The adaptation algorithms are presented in details in the following papers.
  1) BellLab: P. R. Alface, J. Macq and N. Verzijp, "Interactive omnidirectional video delivery: A bandwidth-effective approach," in Bell Labs Technical Journal, vol. 16, no. 4, pp. 135-147, March 2012.
  2) France: X. Corbillon, G. Simon, A. Devlic and J. Chakareski, "Viewport-adaptive navigable 360-degree video delivery," 2017 IEEE International Conference on Communications (ICC), Paris, 2017, pp. 1-7.
  3) ROI: The lowest version is simply chosen for the invisible tiles and the highest possible version for the visible tiles.
  4) ******

## How to use
### Metadata (e.g., Metadata/BR_PSNR_6f1x1_low_delay_300Fr.txt)
* Provide the information (i.e., bitrate, PSNR, and MSE values) of tiles' versions

### Algorithm BellLab (code: Functions/F_BellLab.m)
* Input
  - Double Fh          : Horizontal Field of View
  - Double Fv          : Vertical Field of View
  - Double vp_W        : Width of viewport
  - Double vp_H        : Height of viewport
  - Double erp_W       : Width of video in the original format (i.e., ERP)
  - Double erp_H       : Height of video in the original format (i.e., ERP)
  - Double phi & theta : Position of viewport 
  - Double No_tile     : Total number of tiles
  - Double No_ver      : Number of quality versions
  - Double Bw          : Constrained Bandwidth
  - Array  LB_tile_W   : Lower horizontal boundaries of Tiles (to determine the position of each tile)
  - Array  LB_tile_H   : Lower vertical boundaries of Tiles (to determine the position of each tile)
  - Array  HB_tile_W   : Higher horizontal boundaries of Tiles (to determine the position of each tile)
  - Array  HB_tile_H   : Higher vertical boundaries of Tiles (to determine the position of each tile)
  - Array  Uti         : Marginal utility values of Tiles' Versions 
  - Array  Cti         : Marginal costs of Tiles' Versions 
  
* Output
  - Array  T           : Selected Versions for Tiles
  - Double current_BW  : Total bitrate of selected versions
  - Array  m_ & n_     : Position on ERP format of points on viewport (i.e., Each point (i,j) on viewport coresponds to the point of (m_(i,j),n(i,j)) on ERP format.)
  - Array  P           : Feature of Tiles (i.e., P(i) = 1 -> Tile i is visible | P(i)=0 -> Tile i is invisible) 
  
* Example: BellLab_allBw_8x8.m shows an example of the use of the function F_BellLab. 

### Algorithm France (code: Functions/F_France_Cube.m)
* Input
  - Double Fh          : Horizontal Field of View
  - Double Fv          : Vertical Field of View
  - Double vp_W        : Width of viewport
  - Double vp_H        : Height of viewport
  - Double face_W      : Width of each face in the original format (i.e., Cubemap)
  - Double face_H      : Height of each face in the original format (i.e., Cubemap)
  - Double phi & theta : Position of viewport 
  - Double No_face     : Number of faces
  - Double tile_hori_num     : Number of horizontal tiles
  - Double tile_ver_num      : Number of vertical tiles
  - Double No_ver      : Number of quality versions
  - Double T_e         : Constrained Bandwidth 
  - Array  BR          : Bitrates of Tiles' Versions 
  - Array  MSE         : MSE of Tiles' Versions 
  - Array  LB_tile_W   : Lower horizontal boundaries of Tiles (to determine the position of each tile)
  - Array  LB_tile_H   : Lower vertical boundaries of Tiles (to determine the position of each tile)
  - Array  HB_tile_W   : Higher horizontal boundaries of Tiles (to determine the position of each tile)
  - Array  HB_tile_H   : Higher vertical boundaries of Tiles (to determine the position of each tile)
  
* Output
  - Array  N_ft        : Number of pixels on viewport rounded from each tile (i.e., N_ft(f,t): number of pixels on viewport rounded from tile t of face f)
  - Array  v_sel       : Selected Versions for Tiles
  - Double TB_sel      : Total bitrate of selected versions
  
* Example: France_Cube_allBw_6face.m shows an example of the use of the function F_France_Cube. 

### Algorithm ROI (code: Functions/F_ROI_BG_Cube.m)
* Input
  - Double Fh          : Horizontal Field of View
  - Double Fv          : Vertical Field of View
  - Double vp_W        : Width of viewport
  - Double vp_H        : Height of viewport
  - Double face_W      : Width of each face in the original format (i.e., Cubemap)
  - Double face_H      : Height of each face in the original format (i.e., Cubemap)
  - Double phi & theta : Position of viewport 
  - Double No_face     : Number of faces
  - Double tile_hori_num     : Number of horizontal tiles
  - Double tile_ver_num      : Number of vertical tiles
  - Double No_ver      : Number of quality versions
  - Double Bw          : Constrained Bandwidth 
  - Array  deltaBR     : Bitrate Differences between Tiles' Consecutive Versions 
  - Array  LB_tile_W   : Lower horizontal boundaries of Tiles (to determine the position of each tile)
  - Array  LB_tile_H   : Lower vertical boundaries of Tiles (to determine the position of each tile)
  - Array  HB_tile_W   : Higher horizontal boundaries of Tiles (to determine the position of each tile)
  - Array  HB_tile_H   : Higher vertical boundaries of Tiles (to determine the position of each tile)
  
* Output
  - Array  s_t_        : Number of pixels on viewport rounded from each tile (i.e., N_ft(f,t): number of pixels on viewport rounded from tile t of face f)
  - Array  T           : Selected Versions for Tiles
  - Double sel_BR      : Total bitrate of selected versions	

* Example: ROI_BG_Cube_allBw_6face1x1.m shows an example of the use of the function F_ROI_BG_Cube. 


## Authors

* **Tran Huyen** - *The University of Aizu, Japan* - tranhuyen1191@gmail.com

## Acknowledgments

If you use this source code in your research, please cite

1. The link to this repository.
2. P. R. Alface, J. Macq and N. Verzijp, "Interactive omnidirectional video delivery: A bandwidth-effective approach," in Bell Labs Technical Journal, vol. 16, no. 4, pp. 135-147, March 2012.
3. X. Corbillon, G. Simon, A. Devlic and J. Chakareski, "Viewport-adaptive navigable 360-degree video delivery," 2017 IEEE International Conference on Communications (ICC), Paris, 2017, pp. 1-7.

## License

The source code is only used for non-commercial research purposes.
* If you have any questions, suggestions or corrections, please email to tranhuyen1191@gmail.com. 
