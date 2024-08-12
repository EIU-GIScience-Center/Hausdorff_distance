# -*- coding: utf-8 -*-
"""
Created on Thu Aug  8 09:05:25 2024

@author: Administrator
"""

import utils_data

feat = utils_data.get_sample_feature("cannonball")
print(f"vertices: {len(feat)}")

# # Ramer/Douglas Peucker (1972-3)
from linesimplify import douglas_peucker as dp
errors, sorted_errors = dp.get_errors_sortedErrors(feat) # get displacement distance for each point
# simplify to 500 meet her tolerance
simp_feat = dp.simplify_by_distance_tolerance(feat, errors, 500) # simplify to distance displacement tolerance
print(f"simplified vertices: {len(simp_feat)}")
out_folder = r"sample_data\simplified_dp_500m"
out_file = f"{out_folder}\\cannonball_simplified_dp_500m.shp"
utils_data.create_polyline_shapefile([[simp_feat]], out_file)
# simplify to 10 pct of original vertices
n = int(len(feat)/10)
simp_feat= dp.simplify_by_numPts(feat, n, errors, sorted_errors) # simplify to 10 pts
print(f"simplified vertices: {len(simp_feat)}")
out_folder = r"sample_data\simplified_dp_10pct"
out_file = f"{out_folder}\\cannonball_simplified_dp_10pct.shp"
utils_data.create_polyline_shapefile([[simp_feat]], out_file)


# AREA-PRESERVING SEGMENT COLLAPSE (APSC)
import linesimplify
from linesimplify import apsc
# simp_table=linesimplify.apsc.simplificationTable(pts) # create the simplification table 
# simplifiedline = apsc.simplifiedLine(simp_table,min_pts=10) # simplify to 10 points
# simplifiedline = apsc.simplifiedLine(simp_table,max_disp=50) # simplify to areal displacement tolerance


# # OTHER ALGORITHMS
# # Visvalingam & Whyatt (1993)
# from linesimplify import visvalingam as v
# import numpy as np
# simplifier=v.VWSimplifier(np.array(pts)) # create simplifier object
# simplified_line = simplifier.from_number(10).tolist() # simplify to 10 pts
# simplified_line = simplifier.from_threshold(50).tolist() # simplify to areal displacement tolerance

