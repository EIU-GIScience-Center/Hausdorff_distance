# -*- coding: utf-8 -*-
"""
Created on Sun Jan 28 19:54:22 2024

@author: Administrator
"""

from linesimplify import apsc
from linesimplify import douglas_peucker as dp
import shapefile
import time
import matplotlib.pyplot as plt
import utils_hausdorff
import utils_geom
import polyline_hausdorff
from utils_data import get_sample_feature,copy2clip

def time_test_equal_vertex_counts(feat_a_func,feat_b_func,vertex_counts,brute_force = False):
    msg = f"n\telapsed\tavg_dist"
    print(msg)
    data = msg
    for n in vertex_counts:    
        # Create two versions of a simplified line
        feat_a = feat_a_func(n)
        feat_b = feat_b_func(n)
        # calculate average polyline distance
        start = time.time() 
        avgD = polyline_hausdorff.Hausdorff_average_distance(feat_a,feat_b,brute_force)[1]
        finish = time.time()
        elapsed = finish-start
        msg = f"{n} \t {elapsed} \t {avgD}"
        data = data + "\n" + msg
        print(msg)
    copy2clip(data)
    print("Finished")

def time_test_unequal_vertex_counts(
        feat_a_func,
        feat_b_func,
        vertex_counts,
        brute_force = False):
    msg = "n-m\t" + "\t".join([str(v) for v in vertex_counts])
    print(msg)
    for n in vertex_counts:    
        times = []
        print(f"n: {n}")
        for m in vertex_counts:
            print(f"...m: {m}")
            # Create two versions of a simplified line
            feat_a = feat_a_func(n)
            feat_b = feat_b_func(m)
            # calculate average polyline distance
            start = time.time() 
            avgD = polyline_hausdorff.Hausdorff_average_distance(feat_a,feat_b,brute_force)[1]
            finish = time.time()
            elapsed = finish-start
            times.append(elapsed)
        msg += f"\n{n}\t" + "\t".join([str(t) for t in times])
    copy2clip(msg)
    print("Finished")


def zigzag(n,vertical=False):
    pts = [(i/(n-1),i%2) for i in range(n)]
    if vertical:
        pts = [(y,x) for x,y in pts]
    return pts

def test_same_polyline(feat):
    """
    See if it's really giving a nonzero distance for the same polylines compared with each other'
    """
    simp_table=apsc.simplificationTable(feat) # create the simplification table 
    vertex_counts = [1000,100]
    def get_first_feature(n):
        return apsc.simplifiedLine(simp_table,min_pts=n) 
    time_test_equal_vertex_counts(get_first_feature,get_first_feature,vertex_counts)

def time_test_equal(feat,vertex_counts):
    """ Tests how long it takes to compute average distance between 
        two polylines with equal numbers of vertices """
    # precalculate
    simp_table=apsc.simplificationTable(feat) # create the simplification table 
    errors, sorted_errors = dp.get_errors_sortedErrors(feat) # get displacement distance for each point
    def get_first_feature(n):
        return apsc.simplifiedLine(simp_table,min_pts=n) 
    def get_second_feature(n):
        return dp.simplify_by_numPts(feat, n, errors, sorted_errors) # simplify to 10 pts,[]
    time_test_equal_vertex_counts(get_first_feature,get_second_feature,vertex_counts)

def time_test_unequal(feat,vertex_counts,brute_force = False):
    """ Tests how long it takes to compute average distance between 
        two polylines with equal numbers of vertices """
    # precalculate
    simp_table=apsc.simplificationTable(feat) # create the simplification table 
    errors, sorted_errors = dp.get_errors_sortedErrors(feat) # get displacement distance for each point
    def get_first_feature(n):
        return apsc.simplifiedLine(simp_table,min_pts=n) 
    def get_second_feature(n):
        return dp.simplify_by_numPts(feat, n, errors, sorted_errors) # simplify to 10 pts,[]
    time_test_unequal_vertex_counts(
        get_first_feature,
        get_second_feature,
        vertex_counts,
        brute_force = brute_force)


def time_test_zigzag(vertex_counts):
    def zigzag_vertical(n):
        return zigzag(n,False)
    def zigzag_horizontal(n):
        return zigzag(n,True)
    time_test_equal_vertex_counts(zigzag_vertical,zigzag_horizontal,vertex_counts)

def compare_with_approx(feat,spacings,vertex_ratio = 0.1):
    """
    compare time and accuracy of exact and approximate methods
    """
    msg = "method\tspacing\tnA\tnB\tmax_dist\tfromComp\ttoComp\tavgD\telapsed_sec"
    def add_msg_line(msg,method,spacing,nA,nB,max_dist,fromComp,toComp,avgD,elapsed_sec):
        val = [method,spacing,nA,nB,max_dist,fromComp,toComp,avgD,elapsed_sec]
        msg_line = '\t'.join([f"{x:.6f}" if type(x) == float else str(x) for x in val])
        msg = f"{msg}\n{msg_line}"        
        return msg
    # create simplified feature
    A = feat
    simp_table=apsc.simplificationTable(A,report_interval = 20000) # create the simplification table 
    errors, sorted_errors = dp.get_errors_sortedErrors(A) # get displacement distance for each point
    n = int(len(A) * vertex_ratio)
    B = apsc.simplifiedLine(simp_table,min_pts=n) 

    # # debugging
    # print("debugging...")
    # print("...A > B...")
    # x = polyline_hausdorff.directional_polyline_hausdorff(A,B)
    # print(f'old function A > B: {x}')
    # print("...B > A...")
    # x = polyline_hausdorff.directional_polyline_hausdorff(B,A)
    # print(f'old function B > A: {x}')


    # run exact method
    print("running exact method...")
    start = time.time() 
    max_dist,avgD,fromSeg = polyline_hausdorff.Hausdorff_average_distance(A,B)
    finish = time.time()
    elapsed_sec = finish-start
    msg = add_msg_line(msg,"exact","n/a",len(A),len(B),max_dist,f"s{fromSeg}","n/a",avgD,elapsed_sec)

    # run approximate method
    t = "   ".join([f"{x:.3f}" for x in spacings])
    print(f"spacings to run: {t}")
    for spacing in spacings:
        print(f"running approx method (spacing={spacing})")
        start = time.time() 
        dense_poly = utils_geom.densify(A, spacing)
        max_dist,avgD,fromVertex,nearSeg = approx_polyline_distances(dense_poly,B)
        finish = time.time()
        elapsed_sec = finish-start
        msg = add_msg_line(msg,"approx",spacing,len(dense_poly),len(B),max_dist,fromVertex,nearSeg,avgD,elapsed_sec)
        copy2clip(msg)

        # # debugging
        # A2 = dense_poly[fromVertex-1:fromVertex+3]
        # max_dist,avgD,fromSeg = polyline_hausdorff.Hausdorff_average_distance(A2,B)
        # msg = add_msg_line(msg,"exact","n/a",len(A2),len(B),max_dist,f"s{fromSeg}","n/a",avgD,elapsed_sec)
        # copy2clip(msg)
        
    print(msg)    

def create_example_plots(orig_feat):
    """
    create sample plots for paper demonstrating scenarios with
    differences in computational complexity
    """
    c1 = "blue"
    c2 = "darkorange"
    fig, axes = plt.subplots(1,3,figsize=(7,2))
    plt.axis('equal')
    for ax in axes:
        ax.axis("off")
        plt.axis('equal')
    # show two polylines with equal number of vertices
    n = 25
    simp_table=apsc.simplificationTable(orig_feat) # create the simplification table 
    errors, sorted_errors = dp.get_errors_sortedErrors(orig_feat) # get displacement distance for each point
    feat_a = apsc.simplifiedLine(simp_table,min_pts=n) 
    feat_b = dp.simplify_by_numPts(orig_feat, n, errors, sorted_errors)
    ax = axes[0]
    x,y = zip(*feat_a)
    ax.plot(x,y,color=c1,linewidth=0.5)
    x,y = zip(*feat_b)
    ax.plot(x,y,color=c2,linewidth=0.5)
    ax.set_title('best case:\n equal # vertices', y = -0.3)
    # show one polyline with 10x vertices as the other
    n = 10
    feat_a = apsc.simplifiedLine(simp_table,min_pts=n*10) 
    feat_b = dp.simplify_by_numPts(orig_feat, n, errors, sorted_errors)
    ax = axes[1]
    x,y = zip(*feat_a)
    ax.plot(x,y,color=c1,linewidth=0.5)
    x,y = zip(*feat_b)
    ax.plot(x,y,color=c2,linewidth=0.5)
    ax.set_title('typical case:\n unequal # vertices', y = -0.3)
    # show worst case scenario - 2 zigzags in opposite directions    
    ax = axes[2]
    feat_a = zigzag(21)
    feat_b = zigzag(21,True)
    x,y = zip(*feat_b)
    ax.plot(x,y,color=c1,linewidth=0.5)
    x,y = zip(*feat_a)
    ax.plot(x,y,color=c2,linewidth=0.5)
    ax.set_title('worst case:\n orthogonal zigzags',y = -0.3)
    # display and save    
    
    
    folder = r"C:\CaGIS Board Dropbox\cantaloupe bob\Research\Projects\line difference metrics\Hausdorff\images"
    filename = "features_for_time_testing.png"
    plt.savefig(f"{folder}\\{filename}",dpi=300)
    plt.show()

def approx_polyline_distances(A,B):
    """
    Calculates the average and maximum distances between two polylines

    Parameters
    ----------
    A : List of (x,y) tuples
        The first polyline which we are measuring distances from
    B : List of (x,y) tuples
        The second polyline which we are measuring distances to

    Returns
    -------
    (max,avg) : the two distance metrics

    """
    # create index of segments of second polyline
    idx = utils_hausdorff.seg_idx(B)
    # initialize results
    max_dist = 0
    sum_dist = 0
    # debugging
    nearSeg = -1
    fromVertex = -1
    # loop through vertices of first polyline
    for v in range(len(A)):
        seg = utils_hausdorff.nearSegment(A,B,v,idx)
        comp = (True,seg)
        
        loc = utils_hausdorff.nearLoc(A[v],B,comp)
        d = utils_geom.distance(A[v],loc)
        # print(f"v{v} -> s{seg} ({d} {loc})")
        sum_dist += d
        if d > max_dist:
            max_dist = d
            fromVertex = v
            nearSeg = seg
    avgD = sum_dist / len(A)
    return max_dist,avgD,fromVertex,nearSeg

def accuracy_check():
    # check exact method against approximate method with fine spacing
    # to check for possible errors
    # note that in previous rounds errors were in the range of 100+m
    feat = get_sample_feature("india_bangladesh")
    # create simplified feature
    vertex_ratio = 0.1
    simp_table=apsc.simplificationTable(feat) # create the simplification table 
    errors, sorted_errors = dp.get_errors_sortedErrors(feat) # get displacement distance for each point
    n = int(len(feat) * vertex_ratio)
    B = apsc.simplifiedLine(simp_table,min_pts=n) 

    # loop through segments
    spacing = 2
    tol = 0.1
    report_interval = 1000
    print("POSSIBLE ERRORS:")
    error_count = 0
    messages = []
    for seg in range(len(feat)-1):
        A = feat[seg:seg + 2]
        dense_poly = utils_geom.densify(A,spacing)
        # try to get region of about 1000 vertices around A
        pos = int(n*seg/len(feat))
        B2 = B[max(0,pos - 500):min(len(B)-1,pos + 500)]
        exact_max,exact_avgD,fromSeg = polyline_hausdorff.Hausdorff_average_distance(A,B2)
        approx_max,avgD,fromVertex,nearSeg = approx_polyline_distances(dense_poly,B2)
        # exact max should not be > sum of approx_max,spacing        
        if exact_max > approx_max + spacing + tol or exact_max < approx_max - spacing - tol:
            messages.append(f"seg\t{seg}\tmax\t{exact_max}\tapprox_max\t{approx_max}")
            error_count += 1
            print(f"{error_count} possible errors found...")
            copy2clip("\n".join(messages))
        if seg % report_interval == 0:
            print(f"{seg} segments processed")
    print("\n".join(messages))
    print(f"{error_count} possible errors found")
        
def debug_exact_method():
    ## DEBUGGING
    # create simplified feature
    feat = get_sample_feature("india_bangladesh")
    # create simplified feature
    vertex_ratio = 0.1
    simp_table=apsc.simplificationTable(feat) # create the simplification table 
    errors, sorted_errors = dp.get_errors_sortedErrors(feat) # get displacement distance for each point
    n = int(len(feat) * vertex_ratio)
    B = apsc.simplifiedLine(simp_table,min_pts=n) 
    
    # isolate exact segment from which err occurs
    # compare_with_approx(A,[1000,10])
    seg = 31538
    spacing = 2
    A = feat[seg:seg + 2]
    dense_poly = utils_geom.densify(A,spacing)
    # try to get region of about 1000 vertices around A
    pos = int(n*seg/len(feat))
    B = B[max(0,pos - 500):min(len(B)-1,pos + 500)]
    
    # isolate location of error
    dense_poly = utils_geom.densify(A,1)
    approx_max_dist,approx_avgD,fromVertex,nearSeg = approx_polyline_distances(dense_poly,B)
    max_dist,avgD,fromSeg = polyline_hausdorff.Hausdorff_average_distance(
        A,B,verbose = False)
    print(f'APPROX max_dist: {approx_max_dist} (v{fromVertex} -> s{nearSeg})')
    print(f'EXACT max_dist: {max_dist}')

    messages = []
    messages.append(f"A1\t{A[0][0]}\t{A[0][1]}")
    messages.append(f"A2\t{A[1][0]}\t{A[1][1]}")
    messages.append(f"B997\t{B[997][0]}\t{B[997][1]}")
    messages.append(f"B998\t{B[998][0]}\t{B[998][1]}")
    messages.append(f"B931\t{B[931][0]}\t{B[931][1]}")
    copy2clip("\n".join(messages))
    
    # # visualize situation
    # import matplotlib.pyplot as plt
    # fig, ax = plt.subplots()
    # # plot A
    # x,y = zip(*A)
    # ax.plot(x,y,color = "orange",linewidth = 2,zorder = 1)
    # # x,y = zip(*[B[44]])
    # # ax.scatter(x,y,color = "black", linestyle = "None",zorder = 2)
    # x,y = zip(*B[44:48])
    # ax.plot(x,y,color = "blue",linewidth = 0.5)
    # ax.axis('equal')
    # get exact location of near component
    a = 0
    B_seg_idx = utils_hausdorff.seg_idx(B)
    b = utils_hausdorff.nearSegment(A, B, a, B_seg_idx)
    b2 = utils_hausdorff.nearSegment(A, B, a + 1, B_seg_idx)
    comp, nearloc, d = utils_hausdorff.nearComponent(A, B, a, b)
    comp2, nearloc2, d2 = utils_hausdorff.nearComponent(A, B, a+1, b2)
    comp_list = utils_hausdorff.candidateComponents(
        A,B,a,
        nearloc,nearloc2,
        d,d2,
        B_seg_idx,False
        )
    print(f"B: {len(B)}")
    near_comps = polyline_hausdorff.segment_traversal(A, B, a, comp,d,comp2,d2, comp_list,False)

    comp = max(near_comps,key = lambda x: x[0])
    print(f'MAX COMP: {comp}')
    print("NEAR COMPONENTS:")
    [print(x[:3]) for x in near_comps]
    k=comp[1]
    loc = utils_geom.location(A,0,k)
    
    # comp = comp[2]
    # print(f'comp: {comp}')
    # if comp[0]:
    #     nearloc = B[comp[1]]
    # else:
    #     nearloc = utils_geom.project_pt_to_line(loc, B[comp[1]], B[comp[1]+1])
    # x,y = zip(*[loc])
    # ax.scatter(x,y,color = "black", linestyle = "None",zorder = 2,s = 5)
    # x,y = zip(*[nearloc])
    # ax.scatter(x,y,color = "black", linestyle = "None",zorder = 2,s = 5)
    # x,y = zip(*[loc,nearloc])
    # ax.plot(x,y,color = "red",linewidth = 0.5)
    # filename = r"C:\CaGIS Board Dropbox\cantaloupe bob\Barry\Research\Projects\polyline difference metrics\Hausdorff\images\exact_distance_longer.png"
    # plt.savefig(filename,dpi = 600)
    # plt.show()
    
    # import matplotlib.pyplot as plt
    # fig, ax = plt.subplots()
    # ax.axis('equal')
    # x,y = zip(*B)
    # ax.plot(x,y,color = "black")
    # x,y = zip(*A)
    # ax.plot(x,y,color = "red")
    # x,y = zip(*A[380:382])
    # ax.plot(x,y,color = "turquoise")
    
    # # debugging
    # # create simplified feature
    # feat = feat[80370:80800]
    # vertex_ratio = 0.1
    # A = feat
    # simp_table=apsc.simplificationTable(A,report_interval = 20000) # create the simplification table 
    # errors, sorted_errors = dp.get_errors_sortedErrors(A) # get displacement distance for each point
    # n = int(len(A) * vertex_ratio)
    # B = apsc.simplifiedLine(simp_table,min_pts=n) 

feat = [(0,0),(10,10),(10,0),(20,0)]
feat = get_sample_feature("cannonball")
# print(f"Vertices in input feature: {len(feat)}")
## COMPARISON WITH APPROXIMATE METHODS
y = 10 ** (1/4)
spacings = [1000/(y ** x) for x in range(13)][:1]
compare_with_approx(feat,spacings)
# accuracy_check()
# debug_exact_method()