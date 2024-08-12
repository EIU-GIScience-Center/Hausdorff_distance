# -*- coding: utf-8 -*-
# Usage instructions for polyline Hausdorff and average distance
"""
***********************************************
                 BASIC USAGE
***********************************************
"""

# Import the main module
import polyline_hausdorff as ph

# Input polylines should be lists of (x,y) tuples
A = [(1, 18), (13, 26), (28, 11), (15, 1)]
B = [(1, 17), (10, 23), (10, 19), (13, 19), (13, -2), (26, 8)]

# Compute the hausdorff and average distances. 
d_haus, d_avg, transitions = ph.hausdorff_average_distance(A,B)
print(f"The true Polyline Hausdorff distance is {d_haus:.3f}.")
print(f"The average distance from A to B {d_avg:.3f}.")

"""
***********************************************
   HAUSDORFF LOCATION AND DISTANCE FUNCTION
***********************************************
"""
# 'transitions' contains info on transition points along A
# where nearest component on b changes. It can be used to
# find Hausdorff location and reconstruct distance function

import utils_hausdorff as hu
import utils_geom as gu
from itertools import accumulate

# HAUSDORFF LOCATION
index_max = max(range(len(transitions)), key = lambda x:transitions[x][1])
pos,d,comp = transitions[index_max] # position on A and nearest component on B
a = int(pos) # index of segment on polyline A
k = pos - a # relative position along segment
# get text description of source location
if pos == int(pos):
    src = f"vertex {a} of polyline A"
else:
    src = f"a point {100*k:.1f}% of the way between vertices {a} and {a + 1} of polyline A"
# get text description of target location
if index_max == 0:
    trg = f"{hu.component_label(comp)} of polyline B"
else:
    pos2,d2,comp2 = transitions[index_max - 1]
    trg = f"{hu.component_label(comp)} and {hu.component_label(comp2)} of polyline B"
print(f"The hausdorff distance is from {src} to {trg}.") 

# DISTANCE FUNCTION 
seg_lens = [gu.distance(A[a],A[a+1]) for a in range(len(A)-1)] + [0]
cumulative_lens = [0] + list(accumulate(seg_lens))
dist_func = []
for i in range(len(transitions) - 1):
    pos,d,comp = transitions[i]
    next_pos = transitions[i+1][0]
    a = int(pos)
    while pos < next_pos:
        k = pos - a
        l = cumulative_lens[a] + k * seg_lens[a]
        if k != 0:
            rep = hu.distanceRepresentation(A, B, a, comp)
            d = hu.componentDistance(rep, k, seg_lens[a])
        dist_func.append((l,d))
        pos += 0.05
# dist_func contains tuples of (position along A, distance from B)

"""
***********************************************
    EXAMPLE PLOT OF DISTANCE FUNCTION
***********************************************
"""
import matplotlib.pyplot as plt
# set up figure with 2 plots
fig,axs = plt.subplots(2,figsize = (8,9),height_ratios = [2,1])

# EXTRACT COORDINATES FOR POLYLINE VERTICES AND TRANSITION POINTS
# obtain different sets of vertices
x,y = zip(*A) # polyline A vertices
tx,ty = zip(*[gu.location(A,int(pos),pos - int(pos)) 
              for pos,d,comp in transitions[:-1] 
              if pos != int(pos)])
bx,by = zip(*B)

# EXTRACT POINTS ON DISTANCE FUNCTION FOR VERTICES AND TRANSITION POINTS
v_df = []
t_df = []
for i in range(len(transitions)):
    pos,d,comp = transitions[i]
    a = int(pos)     
    k = pos - a
    if k == 0:
        v_df.append((cumulative_lens[a],d))
    else:
        t_df.append((cumulative_lens[a] + k * seg_lens[a],d))
# separate x & y components for matplotlib
fl,fd = zip(*dist_func)
vl,vd = zip(*v_df)
tl,td = zip(*t_df)

# GET LABEL POSITIONS FOR NEAR COMPONENTS
start_x = 0
df_lbl_x = []
df_comp = [transitions[0][2]]
for t in transitions:
    if t[2] != df_comp[-1]:
        pos,d,comp = t
        a = int(pos)
        k = pos - a
        l = cumulative_lens[a] + k * seg_lens[a]
        df_lbl_x.append((start_x + l)/2)
        start_x = l
        df_comp.append(comp)
df_lbl_x.append((cumulative_lens[-1] + start_x)/2)
df_lbl = [hu.component_label(comp,True) for comp in df_comp]

# PLOT POLYLINES ON FIRST PLOT
ax = axs[0]
# polylines
ax.plot(x,y,
        color = "darkorange",
        label = "A")
ax.plot(x[0],y[0],
        marker = "o",
        markersize=7,
        color = "black",
        linestyle = "None",
        label = "start")
ax.plot(x[-1],y[-1],
        marker = "X",
        markersize=7,
        color = "black",
        linestyle = "None",
        label = "finish")
ax.plot(x[1:-1],y[1:-1],
        linestyle = "None",
        marker = "o",
        markersize = 5,
        markerfacecolor = "darkorange",
        markeredgecolor = "black",
        label = "vertex")
ax.plot(tx,ty,
        linestyle = "None",
        marker = "o",
        markersize = 3,
        markerfacecolor = "black",
        markeredgecolor = "black",
        label = "transition")
ax.plot(bx,by,"dodgerblue",label = "B")
ax.plot(bx,by,
        linestyle = "None",
        marker = "o",
        markersize = 12,
        markerfacecolor = "white",
        markeredgecolor = "dodgerblue",
        label = "vertex")

for curX,curY,lbl in zip(bx,by,range(len(bx))):
    ax.text(curX,curY,lbl,ha = "center",va = "center_baseline")
# refine plot
ax.set_title("Input Polylines")
ax.axis("equal")
ax.legend()

# PLOT DISTANCE FUNCTION ON SECOND PLOT
# start plotting
ax = axs[1]
ax.plot(fl,fd,color = "gray")
ax.plot(tl,td,
        linestyle = "None",
        marker = "o",
        markersize = 3,
        markerfacecolor = "black",
        markeredgecolor = "black",
        label = "transition points")
ax.plot(vl[1:-1],vd[1:-1],
        linestyle = "None",
        marker = "o",
        markersize = 5,
        markerfacecolor = "darkorange",
        markeredgecolor = "black",
        label = "vertices")
ax.plot(fl[0],fd[0],
        marker = "o",
        markersize=7,
        color = "black",
        linestyle = "None",
        label = "start")
ax.plot(fl[-1],fd[-1],
        marker = "X",
        markersize=7,
        color = "black",
        linestyle = "None",
        label = "finish")
# Labels
# Label nearest component on x-axis
ax.set_xticks(tl,[])
ax.set_xticks(df_lbl_x,df_lbl,minor = True)
ax.tick_params(axis='x', which="minor",length=0)
ax.set_xlabel("position along A \nlabels indicate nearest component on B")
# Y-axis and title
ax.set_ylabel("distance to B")
ax.set_title("Distance Function (from A to B)")
plt.show()
    
# import math

# # plot the original polylines
# x,y = zip(*A)
# plt.plot(x,y,linewidth=3,c='darkorange',marker='o',markersize=5,mfc='white',mec='gray')
# x,y = zip(*B)
# plt.plot(x,y,linewidth=3,c='blue',marker='o',markersize=5,mfc='white',mec='gray')
# plt.axis('equal')

# # plot the source location
# srcloc = hausdorff_info[1]  
# plt.plot(srcloc[0],srcloc[1],marker='X',mec='red',mfc='red',markersize=12)

# # plotting text and arrows is complicated - let's make a function for it
# def showDist(srcloc,trgloc,df_lbl,augment=0):
#     # draw arrow from source to target
#     plt.arrow(srcloc[0],srcloc[1],(trgloc[0]-srcloc[0]),(trgloc[1]-srcloc[1]),color='black',head_width=0.25,head_length=0.25,overhang=0.15,length_includes_head = True)
#     # calculate rotation parameters for text
#     r = math.degrees(math.atan2(trgloc[1]-srcloc[1],trgloc[0]-srcloc[0]))% 360
#     if 100 < r < 260:
#         r -= 180
#     # calculate text placement
#     offset = 0.5
#     textx,texty = -1*offset*math.sin(math.radians(r)) + (srcloc[0]+trgloc[0])/2, offset*math.cos(math.radians(r)) + (srcloc[1]+trgloc[1])/2    
#     # plot text
#     plt.text(textx,texty,df_lbl,ha='center',va='center',fontsize=15+augment,rotation=r)
#     # calculate distance
#     dx2 = (trgloc[0] - srcloc[0])**2 
#     dy2 = (trgloc[1] - srcloc[1])**2 
#     d = math.sqrt(dx2+dy2)
#     # show on plot    offset=-0.5
#     textx,texty = -1*offset*math.sin(math.radians(r)) + (srcloc[0]+trgloc[0])/2, offset*math.cos(math.radians(r)) + (srcloc[1]+trgloc[1])/2    
#     plt.text(textx,texty,"{:.2f}".format(d),ha='center',va='center',fontsize=12+augment,rotation=r)    

# # plot the target locations, with text and arrows
# trglocs = hausdorff_info[2]
# for trgloc in trglocs:
#     showDist(srcloc,trgloc,"polyline hausdorff")




# """
# ***********************************************
# COMPARE WITH SCIPY
# ***********************************************
# """
# try:
#     from scipy.spatial.distance import directed_hausdorff
#     import numpy as np
#     # create numpy arrays
#     u = np.array(A)
#     v = np.array(B)
#     # calculate hausdorff distance and get participating vertices
#     d1,a1,b1 = directed_hausdorff(u,v)
#     d2,b2,a2 = directed_hausdorff(v,u)
#     if d1 > d2:
#         d = d1
#         srcloc = A[a1]
#         trgloc = B[b1]
#     else:
#         d = d2
#         srcloc = B[b2]
#         trgloc = A[a2]
#     # print results
#     print("SciPy Hausdorff distance is {:.2f}.".format(d))
#     # show graphically
#     showDist(srcloc,trgloc,"scipy")
# except ImportError:
#     print("Install scipy & numpy to see a comparison with SciPy's Directed Hausdorff distance.")




# """
# ***********************************************
# COMPARE WITH SHAPELY
# ***********************************************
# """
# import utils_hausdorff as uh
# # shapely appears to calculate the largest vertex-to-edge gap
# try:
    
#     from shapely.geometry import LineString
#     # create linestring objects
#     lineA = LineString(A)
#     lineB = LineString(B)
#     # calculate hausdorff distance
#     shapely_hausdorff = lineA.hausdorff_distance(lineB)
#     # print results
#     print("Shapely Hausdorff distance is {:.2f}.".format(shapely_hausdorff))
#     # show graphically
#     srcloc = B[2] 
#     trgloc = uh.nearLoc(srcloc, A, (True,1)) 
#     showDist(srcloc,trgloc,"shapely")
# except ImportError:
#     print("Install shapely to see a comparison with shapely's Hausdorff distance.")



    
# plt.show()
# print("finished.")