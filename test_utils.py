# -*- coding: utf-8 -*-
"""
test data for testing and illustrating Hausdorff distance
"""

import matplotlib.pyplot as plt
import math
import utils_geom as g
import polyline_hausdorff as h
import utils_hausdorff as hu

def case(casenum):
    """
    Creates a test case consisting of two polylines

    Parameters
    ----------
    casenum : int
        Number between 0 and count -1.

    Returns
    -------
    Two polylines, or if casenum == -1 returns the number of available cases

    """
    cases = []
    
    # 0 input lines from Hangouet 1995, fig. 13
    # slightly modified so that solution is as indicated by Hangouet
    A = [(0,0),(0.44,-1.99),(5,0.95)]
    B = [(3.5,-1.2),(3.19,-0.68),(5,0.3),(5.61,1.77),(3.19,0.59),(0.18,0.59),(0.77,-1.3)]
    B.reverse()
    cases.append((A,B))
    
    # 1 input lines from Hangouet 1995, fig. 13
    # modified again to my liking
    A = [(0.44,-1.99),(5,0.95)]
    B = [(3,-1.3),(3.5,-1.2),(3.19,-0.68),(5,0.3),(5.61,1.77),(3.19,0.59),(0.8,0.59),(0.77,-1.3)]
    B.reverse()
    cases.append((A,B))
    
    # 2 parallel lines
    A = [(0,0),(3,3),(12,2),(14,5)]
    B = [(0,12),(4,6),(13,5),(9,6)]
    cases.append((A,B))
    
    # 3 coincident lines
    A = [(0,0),(4,6),(12,2),(14,5),(17,7)]
    B = [(0,10),(4,6),(12,2),(9,7)]
    cases.append((A,B))
    
    # 4 overlapping lines
    A = [(0,0),(4,6),(12,2),(14,5),(17,7)]
    B = [(0,10),(6,5),(10,3),(9,7)]        
    cases.append((A,B))

    # 5 two parallel segments in a row
    A = [(0,0),(4,6),(14,1),(16,3),(19,5)]
    B = [(0,10),(6,6),(8,5),(10,4),(9,7)]        
    cases.append((A,B))
    
    # 6 two overlapping segments in a row
    A = [(0,0),(4,6),(14,1),(16,3),(19,5)]
    B = [(0,10),(6,5),(8,4),(10,3),(9,7)]        
    cases.append((A,B))

    # 7 situation where we have to double back on a component
    A = [(5,7),(25,14)]
    B = [(17,6),(16,9),(20,7),(27,8),(25,17),(6,12)]
    B.reverse()
    cases.append((A,B))    

    # 8 coincident vertices
    A = [(1,0),(10,1)]
    B = [(0,3),(2,2),(5,4),(3,7),(5,9),(7,6),(5,4),(11,3)]
    cases.append((A,B))
    
    # 9 coincident vertices
    A = [(1,0),(10,1)]
    B = [(0,6),(2,5),(5,4),(3,7),(5,9),(7,6),(5,4),(11,6)]
    cases.append((A,B))
    
    # 10 lots of singularities
    A = [(5,3),(25,3)]
    B = [(5,3),(2,3),(2,6),(8,0),(10,0),(12,0),(14,0),(14,8),(12,8),(12,6),(10,10),(16,10),(16,3),(16,0),(18,0),(20,2),(26,2),(25,7),(20,4)]
    cases.append((A,B))
    
    # 11 Hangouët revised - journal article version
    A = [(0, 0), (5, 0)] # plotting code assumes a horizontal line beginning at (0,0)
    B = [(2.4, -1), (2.5, -0.4), (4.5, -0.9), (6.4, 0.8), (4, 0.5), (2, 2.3), (0.5, 0.5)]
    cases.append((A,B))
    
    # * cleanup *
    if casenum == -1:
        return len(cases)
    else:
        return cases[casenum]
    
def plot_hausdorff_solution(
        A, a, B,
        title = None, 
        imagefile = None, 
        show_labels = True,
        B_vertex_size = 10,
        B_show_vertex_labels = True,
        B_show_segment_labels = False,
        A_show_near_labels = True,
        k_use_hausdorff = True,
        k = 0.5,
        k_resolution = 0.005, # spacing of points along segment
        k_buffer = 0.0, # plotted overlap between adjacent near components
        stop_distance_function_at_k = False,
        dpi = 90,
        label_start_finish = False
        ):
    # Plots the input lines with the location (of the far Hausdorff location) 
    # and a circle at distance d from the location

    # turn off interactive mode
    plt.ioff()

    # determine Hausdorff solution
    B_idx = hu.seg_idx(B)
    startNearSeg = hu.nearSegment(A, B,a,B_idx)
    endNearSeg = hu.nearSegment(A, B, a+1,B_idx)
    startcomp, startloc, startd = hu.nearComponent(A, B, a, startNearSeg)
    endcomp, endloc, endd = hu.nearComponent(A, B, a+1, endNearSeg)
    candidates = hu.candidateComponents(A,B,a,startloc,endloc,startd,endd,B_idx)
    near_comps = h.segment_traversal(A, B, a, startcomp, startd, endcomp, endd, candidates)
    d_Hausdorff,k_Hausdorff,hausdorff_comps = h.segment_hausdorff(near_comps)
    len_a = g.distance(A[a],A[a+1])
    
    # calculate average distance
    avgD = h.segment_dist_integral(near_comps, len_a, verbose = True)/len_a
    
    # determine which k-value will be the focus
    if k_use_hausdorff:
        k = k_Hausdorff
        d = d_Hausdorff
    else:
        # find distance at k
        # check k is in unit interval
        k = max(0,k)
        k = min(1,k)
    # get component nearest to k
    # near_comps is list of (d,k,comp,rep)        
    i = 0
    while (i < len(near_comps)-1) and (near_comps[i+1][1] < k):
        i += 1 
    near_comp = near_comps[i][2]
    d = hu.componentDistance(near_comps[i][3], k, len_a)

    # get useful geometry
    ax,ay=list(zip(*A[a:a+2]))
    bx,by=zip(*B)
    adx = ax[1]-ax[0]
    ady = ay[1]-ay[0]
    theta = math.degrees(math.atan2(ady, adx))
        
    d0 = hu.componentDistance(near_comps[0][3], 0, len_a)
    d1 = hu.componentDistance(near_comps[len(near_comps)-1][3], 1, len_a)
    if A[a][0] < A[a+1][0]:
        plot_left = A[a][0] - d0
        plot_right = A[a+1][0] + d1
    else:
        plot_left = A[a+1][0] - d1
        plot_right = A[a][0] + d0
    if A[a][1] < A[a+1][1]:
        plot_bottom = A[a][1] - d0
        plot_top = A[a+1][1] + d1
    else:
        plot_bottom = A[a+1][1] - d1
        plot_top = A[a][1] + d0
    
    min_x = min([v[0] for v in B])
    max_x = max([v[0] for v in B])
    min_y = min([v[1] for v in B])
    max_y = max([v[1] for v in B])
    plot_left = min(plot_left,min_x)
    plot_right = max(plot_right,max_x)
    plot_bottom = min(plot_bottom,min_y)
    plot_top = max(plot_top,max_y)
    buffer_pct = 5
    spanx = (plot_right - plot_left) * buffer_pct/100
    spany = (plot_top - plot_bottom) * buffer_pct/100
    plot_left = plot_left - spanx
    plot_right = plot_right + spanx
    plot_bottom = plot_bottom - spany
    plot_top = plot_top + spany
    
    # set up plot
    plt.figure(dpi = dpi)
    plt.axis([plot_left,plot_right,plot_bottom,plot_top])
    plt.gca().set_aspect('equal', adjustable='box')

    # highlight nearest component
    b = near_comp[1]
    if near_comp[0]: # component is a segment
        x = [B[b][0],B[b+1][0]]
        y = [B[b][1],B[b+1][1]]
        plt.plot(x,y,color = "cornflowerblue", linewidth = 2.6)
    else:
        x = [B[b][0]]
        y = [B[b][1]]
        plt.plot(x,y,marker = 'o', 
                    markersize = B_vertex_size+3, 
                    linewidth = 2, 
                    markeredgecolor = 'navy',
                    markerfacecolor = None)

    # plot input lines
    plt.plot(ax,ay, color = "sandybrown", linewidth = 3,
             solid_capstyle = 'butt')
    plt.plot(bx,by, color = "cornflowerblue", linewidth = 0.6, marker = 'o', 
             markersize = B_vertex_size,
             markerfacecolor = 'white')
    
    # label start and finish of segment a 
    off = 0.025
    x,y = g.location(A,a,0-off)
    plt.text(x,y,'start', 
             horizontalalignment = 'right',
             verticalalignment = 'center',
             rotation = theta, 
             rotation_mode = "anchor")
    x,y = g.location(A,a,1+off)
    plt.text(x,y,'finish', 
             horizontalalignment = 'left',
             verticalalignment = 'center',
             rotation = theta, 
             rotation_mode = "anchor")
    
    # label vertices
    axis = plt.gca()
    if B_show_vertex_labels:
        if B_vertex_size >= 10:
            off = (0,0)
            fs = min(B_vertex_size-3,8)
        else:
            off = (0,0.1)
            fs = 10
        for i in range(len(B)):
            b = B[i]
            axis.text(b[0]+off[0],b[1]+off[1],"{}".format(i), 
                      fontsize = fs, 
                      fontweight = 'bold', 
                      horizontalalignment = 'center',
                      verticalalignment = 'center')
            
                    
    # label segments
    if B_show_segment_labels:
        off = (-0.1,-0.1)
        for i in range(len(B)-1):
            b0 = B[i]
            b1 = B[i+1]
            label_loc = (off[0]+(b0[0]+b1[0])/2,off[1]+(b0[1]+b1[1])/2)
            plt.text(label_loc[0],label_loc[1],"s{}".format(i), fontsize = 7, horizontalalignment = 'center')

    # Mark divisions along a
    p = 0.01 # determines how long the division marks are
    # p = p*dpi/90 # adjust for figure resolution
    mdx = -ady * p 
    mdy = adx * p
    for i in range(len(near_comps)):
        k1 = near_comps[i][1]
        if (stop_distance_function_at_k == False) or (k1 <= k):
            x,y = g.location(A,a,k1)
            plt.plot([x-mdx,x+mdx],[y-mdy,y+mdy],color = "gray", linewidth = 0.8)        
            # label section
            if i < len(near_comps) - 1:
                k2 = near_comps[i+1][1]
                x2,y2 = g.location(A,a,k2)
                x,y = (x+x2)/2,(y+y2)/2
                if A_show_near_labels:
                    label = hu.component_label(near_comps[i][2])
                    if near_comps[i][2][0] == True: # segment
                        fs = 7
                        fw = 'normal'
                    else:
                        fs = 10
                        fw = 'bold'
                    plt.text(x+mdx,y+mdy,label,
                             fontsize = fs,
                             fontweight = fw,
                             rotation = theta,
                             rotation_mode = 'default',
                             horizontalalignment = 'center')
    
    # plot distance function
    len_a = g.distance(A[a],A[a+1])
    rescale_factor = 1/((mdx**2+mdy**2)**0.5)
    dx = mdx * rescale_factor
    dy = mdy * rescale_factor
    if dy < 0:
        dx = -dx
        dy = -dy
    for i in range(len(near_comps)-1):
        # near_comps is list of (d,k,comp,rep)     
        k1 = near_comps[i][1]-k_buffer
        k2 = near_comps[i+1][1]+k_buffer
        x,y = g.location(A,a,k1)
        bcomp =  near_comps[i][2]
        rep = hu.distanceRepresentation(A, B, a, bcomp)
        interval = (k1,k2)
        n = max(5, int((k2-k1)/k_resolution))
        if bcomp[0] == True:
            kds = hu.segDistAcross(rep, interval, len_a, n)
        else:
            kds = hu.vertDistAcross(rep, interval, len_a, n)
        distance_line = []
        for k1,distance in kds:
            if (stop_distance_function_at_k == False) or (k1 <= k):
                x,y = g.location(A,a,k1)
                x += distance * dx
                y += distance * dy
                distance_line.append((x,y))
        if len(distance_line) > 0:
            x,y = zip(*distance_line)
            plt.plot(x,y,'--',linewidth=0.8,c='black',marker='o',markersize=0,mfc='white',mec='gray')
            
    # plot circle around Hausdorff solution point
    degs = list(range(361))
    rads = [math.radians(x) for x in degs]
    x_mark,y_mark = g.location(A,a,k)
    pts = [(x_mark+d*math.sin(r),y_mark+d*math.cos(r)) for r in rads]
    x,y = zip(*pts)
    plt.plot(x,y,color='grey')

    # plot distance radius at same point
    d_rad = (x_mark + d * dx, y_mark + d * dy)
    pts = [(x_mark, y_mark),d_rad]
    x,y = zip(*pts)    
    plt.plot(x,y,color='pink')    
    
    # plot red circle at current k-value
    plt.plot(x_mark,y_mark,'ro')

    # write plot title
    if title == None:    
        title = ''
    else:
        title = title + "          "
    title = title + "A: orAnge     B: Blue\n"
    
    # construct plot title
    # msg = " -> ".join([hu.component_label(near_comps[i][2])for i in range(len(near_comps))])
    if k_use_hausdorff:    
        msg = "Hausdorff dist: {:.3f} (k = {:.3f})".format(d,k)
        msg += "    Avg dist: {:.3f}".format(avgD)
    else:
        msg = "dist: {:.3f} (k = {:.3f})".format(d,k)
    title = title + msg
    
    plt.suptitle(title, fontsize = 9)

    # show plot and/or save to image file
    if imagefile != "" and imagefile != None:
        # print("**********************")
        # print("saving to {}".format(imagefile))
        plt.savefig(imagefile)
        plt.clf()
        plt.close
    else:
        # print("**********************")
        # print("showing in interpreter:")
        plt.show()
        plt.clf()
        plt.close()