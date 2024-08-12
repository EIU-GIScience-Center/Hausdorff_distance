# -*- coding: utf-8 -*-
"""
Calculates the Hausdorff distance between two polylines
following the method described in:
J.F. Hangouet (1995). Computation of the Hausdorff distance between plane 
vector polylines. Proceedings, Auto Carto 12, pp. 1-10.

A joint project of the Spring 2021 semester GEO 4910 class at Eastern Illinois 
University.

Team members:
    Luke Jansen
    Tanner Jones
    Farouk Olaitan
    Megshi Thakur
    Barry Kronenfeld (Instructor)
    
License: MIT License

If you publish or create derivate software using this code, acknowledgement 
would be appreciated.

********************************************************************

CONVENTIONS:
    The two polylines are referred to as A and B.

    Lowercase a or b indicates the index of either a vertex or segment on 
    A or B.

    For a polyline with n vertices, vertices are indexed 0 to n-1 and segments 
    are indexed 0 to n-2

DEFINITIONS:
    k-value:        float
        The "relative location" along a segment of polyline A.
        - At the start of the segment, k=0
        - At the end of the segment, k=1
    q-value:        float
        The "relative distance" to a segment of polyline A, i.e. the 
        perpendicular distance to the segment divided by the length of the 
        segment.
    component:      (bool, int)
        A tuple of two values identifying a component of polyline B, consisting of:
            bool:   If True the component is a segment, if False it is a vertex
            int:      The index of the vertex or segment on polyline B
    distance representation:      (bool, float, float)
        A tuple of three values containing information needed to compute the 
        distance from any point on a segment "a" of A to component b of B:
            bool:   If True b is a segment, if False it is a vertex
            float:  The k-value of the projection of vertex b onto segment a, or
                    The k-value of the point of intersection of segment b with segment a, or
                    None if b is a segment parallel to a
            float:  The sin of the angle between a & b (if b is a non-parallel segment)
                    The q-distance from b to a (if b is a vertex or segment parallel to a)
"""

import utils_geom as g
import utils_hausdorff as hu
import integrals

def polyline_hausdorff(A,B):
    """
    ***deprecated, should use Hausdorff_average_distance function***
    
    Computes the Hausdorff distance between two polylines, i.e. the maximum
    distance from any point on either polyline to the nearest point on the
    other polyline.

    Parameters
    ----------
    A : [(x,y),...] list of tuples 
        The coordinates of a polyline.
    B : [(x,y),...] list of tuples 
        The coordinates of another polyline.
    
    Returns
    ----------
    h_dist : float
        The Hausdorff distance between the polylines.
    srcloc : (float,float)
        The location of the source of the Hausdorff distance.
        The source location is the location on either polyline that is furthest 
        from the nearest location on the other polyline.
    srcline : int
        Flag to indicate the source polyline. 
        If 0, the source location is on polyline A. 
        If 1, the source location is on polyline B.
    srccomp : (bool,int)
        The index of the component (vertex or segment) of the source polyline 
        containing the source location. 
        If bool is True then srccomp is a segment, otherwise it is a vertex.
    trglocs : [(x,y),xy]
        The location of the target of the Hausdorff distance.
    trgcomps : [(bool, int)]
        List of components of the target polyline that is equally close to 
        the source location. If the source component is a segment then then 
        there will usually be two equidistant nearest locations on the target 
        polyline.    
    """
    # check in both directions
    h_dist_a,srcloca,srccompa,trgcompsa = directional_polyline_hausdorff(A,B)
    h_dist_b,srclocb,srccompb,trgcompsb = directional_polyline_hausdorff(B,A)
    # determine which is source and which is target    
    if h_dist_a > h_dist_b:
        srcloc = srcloca
        srcline = "A"
        trgline = B
        srccomp = srccompa
        trgcomps = trgcompsa
        h = h_dist_a        
    else:
        srcloc = srclocb
        srcline = "B"
        trgline = A
        srccomp = srccompb
        trgcomps = trgcompsb
        h = h_dist_b
    # determine target location(s)
    trglocs = [hu.nearLoc(srcloc,trgline,trgcomp) for trgcomp in trgcomps]
    
    # report results
    return h,srcloc,srcline,srccomp,trglocs,trgcomps


# from dataclasses import dataclass, field

# @dataclass
# class HausdorfResults:
#     H: float = field(description='')
    
def hausdorff_average_distance(A,B,verbose = False,brute_force = False):
    """
    Computes the one-way Hausdorff and average distance from A to B

    Parameters
    ----------
    A : [(x,y),...] list of tuples 
        The coordinates of the main polyline.
    B : [(x,y),...] list of tuples 
        The coordinates of the other polyline.
    Returns
    ----------
    H : float
        The unidirectional Hausdorff distance between the polylines.
    avgD : float
        The unidirectional average distance between the polylines.
    dist_func : list of (pos,d,comp)
        pos : float
            position along A; int portion is vertex, decimal is k-value
        d : float
            distance to nearest point on B
        comp : (bool, int)
            nearest component on B for following section of distance function
    """

    # initialize component and distance lists
    vNearComp = []
    vNearLoc = []
    vertDist = []
    # create index of B segments
    B_seg_idx = hu.seg_idx(B)
    
    # get distance from each vertex on A to its nearest component on B
    for a in range(len(A)):
        b = hu.nearSegment(A, B, a, B_seg_idx)
        comp, nearloc, d = hu.nearComponent(A, B, a, b)
        vNearComp.append(comp)
        vNearLoc.append(nearloc)
        vertDist.append(d)

    # traverse segments on A tracking nearest components on B
    tot_area = 0
    tot_len = 0
    d_Hausdorff = 0
    # h_from_seg = -1
    dist_func = []
    for a in range(len(A)-1):
        comp_list = hu.candidateComponents(
            A,B,a,
            vNearLoc[a],vNearLoc[a+1],
            vertDist[a],vertDist[a+1],
            B_seg_idx,brute_force
            )
        near_comps = segment_traversal(A, B, a, vNearComp[a],vertDist[a],vNearComp[a+1],vertDist[a+1], comp_list,verbose)
        dist_func += [(a+k,d,comp) for d,k,comp,rep in near_comps]
        # update maximum (Hausdorff) distance
        seg_max = max(x[0] for x in near_comps)
        if seg_max > d_Hausdorff:
            d_Hausdorff = seg_max
            # h_from_seg = a
        # increment total length and total area for average distance computation
        L = g.distance(A[a], A[a+1])
        tot_len += L
        area = segment_dist_integral(near_comps, L)
        tot_area += area

    # handle final vertex
    a = len(A) - 1
    d = vertDist[a]
    k = 0
    comp = vNearComp[a]
    dist_func.append((a+k,d,comp))
    if d > d_Hausdorff:
        d_Hausdorff = d

    d_avg = tot_area / tot_len
    
    return d_Hausdorff, d_avg, dist_func


def near_components(A,B):
    """
    Computes the sequence of nearest components on B from locations on polyline A

    Parameters
    ----------
    A : [(x,y),...] list of tuples 
        The coordinates of the main polyline.
    B : [(x,y),...] list of tuples 
        The coordinates of the other polyline.
    Returns
    ----------
    [a,k,d,comp,rep] : float
        The segment on A, k-value and distance at the start of the section, 
        and nearest component on B and its nearest distance representation
        for successive sections on polyline A
    """

    # initialize component and distance lists
    vNearComp = []
    vNearLoc = []
    vertDist = []
    # create index of B segments
    B_seg_idx = hu.seg_idx(B)

    # get distance from each vertex on A to its nearest component on B
    for a in range(len(A)):
        b = hu.nearSegment(A, B, a, B_seg_idx)
        comp, nearloc, d = hu.nearComponent(A, B, a, b)
        vNearComp.append(comp)
        vNearLoc.append(nearloc)
        vertDist.append(d)

    segsToCheck = [a for a in range(len(A)-1) if hu.checkSegment(vNearComp[a],vNearComp[a+1])]
    # get distance from segments on A to nearest components on B
    a_k_d_comp_rep = []
    for a in segsToCheck:
        comp_list = hu.candidateComponents(
            A,B,a,
            vNearLoc[a],vNearLoc[a+1],
            vertDist[a],vertDist[a+1],
            B_seg_idx,False
            )
        near_comps = segment_traversal(A, B, a, vNearComp[a],vertDist[a],vNearComp[a+1],vertDist[a+1], comp_list)
        # near_comps : list of (d,k,comp,rep)
        for d,k,comp,rep in near_comps:
            a_k_d_comp_rep.append((a,k,d,comp,rep))
    return a_k_d_comp_rep

def directional_polyline_hausdorff(A,B,brute_force = False):
    """
    ***deprecated, we shouldn't need this anymore***
    
    Computes the one-way Hausdorff distance from A to B, i.e. the maximum
    distance from any point on A to the nearest point on B.

    Parameters
    ----------
    A : [(x,y),...] list of tuples 
        The coordinates of the main polyline.
    B : [(x,y),...] list of tuples 
        The coordinates of the other polyline.
    Returns
    ----------
    H : float
        The unidirectional Hausdorff distance between the polylines.
    loc : (float, float)
        The location on A that is the source of the Hausdorff distance.
    srcComp : (bool,int)
        The component of A that is the source of the Hausdorff distance.
    trgComp1 : (bool, int)
        A component of B that is the target of the Hausdorff distance.
    trgComp2 : (bool, int)
        Another component of B that is the target of the Hausdorff distance, or
        None if there is only one target component.
    """

    # initialize component and distance lists
    vNearComp = []
    vNearLoc = []
    vertDist = []
    # create index of B segments
    B_seg_idx = hu.seg_idx(B)

    # get distance from each vertex on A to its nearest component on B
    for a in range(len(A)):
        b = hu.nearSegment(A, B, a, B_seg_idx)
        comp, nearloc, d = hu.nearComponent(A, B, a, b)
        vNearComp.append(comp)
        vNearLoc.append(nearloc)
        vertDist.append(d)

    segsToCheck = [a for a in range(len(A)-1) if hu.checkSegment(vNearComp[a],vNearComp[a+1])]
    # get distance from segments on A to nearest components on B
    d_k_comps = []
    for a in segsToCheck:
        comp_list = hu.candidateComponents(
            A,B,a,
            vNearLoc[a],vNearLoc[a+1],
            vertDist[a],vertDist[a+1],
            B_seg_idx,brute_force
            )
        near_comps = segment_traversal(A, B, a, vNearComp[a],vertDist[a],vNearComp[a+1],vertDist[a+1], comp_list)
        d,k,comps = segment_hausdorff(near_comps)
#         d,k,comps = segment_hausdorff(A, B, a, vNearComp[a],vertDist[a],vNearComp[a+1],vertDist[a+1], comp_list)
        d_k_comps.append((d,k,comps))
    # Hausdorff distance is maximum of distances from vertices and segments
    use_vert = True
    maxVertDist = max(vertDist)
    if len(d_k_comps)!=0:
        maxSegInfo = max(d_k_comps, key = lambda x: x[0])
        if maxSegInfo[0] >= maxVertDist:
            use_vert = False
    if use_vert:
        a = vertDist.index(maxVertDist)
        loc = A[a]
        srcComp = (False,a)
        trgComps = [vNearComp[a]]
        H = maxVertDist
    else: # use segment
        a = segsToCheck[d_k_comps.index(maxSegInfo)]
        d,k,trgComps = maxSegInfo
        loc = g.location(A,a,k)
        srcComp = (True,a)
        H = d
    
    return H,loc,srcComp,trgComps

def segment_hausdorff(seg_traversal, verbose = False):
    """
    Identifies the Hausdorff distance from a traversal of a segment.
    
    Parameters
    ----------
    seg_traversal : [(float,float,comp,dist_rep),...] list of tuples 
        list of (d, k, component, distance representation) along traversal of a
    Returns
    ----------
    d : float
        The directed Hausdorff distance from segment a to polyline B.
    k : float
        The k-value of the source of the Hausdorff distance along segment a.        
    max_comps : [(bool, int)]
        List of one or two components on B that is the target of the Hausdorff distance
    """
    w = 0
    max_dist = 0
    # Loop through traversal to get switch point with maximum distance
    if verbose:
        print("\npositions along A where nearest component on B changes:")
    for i in range(len(seg_traversal)):
        d,k,comp,dist_rep = seg_traversal[i]
        if verbose:
            print("{:.3f}: {} d={:.3f}".format(k,hu.component_label(comp),d))
        if d > max_dist:
            w = i
            max_dist = d
    # get k-value of maximum distance
    k = seg_traversal[w][1]
    # get one component if winner is zero, two components otherwise        
    comps = [seg_traversal[w][2]]
    if w > 0 and seg_traversal[w][2] != seg_traversal[w-1][2]:
        comps.append(seg_traversal[w-1][2])
    if verbose:
        print("red dot is furthest point on A from B")
        print("nearest components on B to red dot: {}".format(",".join([hu.component_label(comp)for comp in comps])))
    return max_dist,k,comps

def segment_dist_integral(seg_traversal,L, verbose = False):
    """
    Compute the integral of the distance function from a traversal of a segment.
    
    Parameters
    ----------
    seg_traversal : [(float,float,comp,dist_rep),...] list of tuples 
        list of (d, k, component, distance representation) along traversal of a
    L : float
        length of traversed segment
    Returns
    ----------
    area : float
        The total area under the distance function from segment a to polyline B.
    """
    # Loop through traversal to get switch point with maximum distance
    if verbose:
        print("\nsections of A nearest to components of B:")
    tot_area = 0
    for i in range(len(seg_traversal)-1):
        d,k0,comp,dr = seg_traversal[i]
        d1,k1,comp1,dr1 = seg_traversal[i+1]
        # use first distance representation
        area = integrals.component_integral(L, k0, k1, dr)        
        tot_area += area
        if verbose:
            print("{:.3f} - {:.3f} {} A={:.3f}".format(k0,k1,hu.component_label(comp),area))
    return tot_area

def segment_traversal(A,B,a,startComp,startDist,endComp,endDist,comp_list,verbose = False,max_iterations = 1000000):
    """
    Computes the furthest distance from any point on segment a of polyline A 
    to the nearest point on polyline B

    Parameters
    ----------
    A : [(x,y),...] list of tuples 
        The coordinates of the main polyline.
    B : [(x,y),...] list of tuples 
        The coordinates of the other polyline.
    a : int
        The index of a segment on polyline A.
    startComp : comp
        Component of B nearest to vertex a.
    startDist : float
        Distance from vertex a to startComp
    endComp : comp
        Component of B nearest to vertex a+1.
    endDist : float
        Distance from vertex a+1 to endComp
    comp_list : [comp]
        List of all candidate components of B for the Hausdorff distance.
    verbose : bool
        If True, will print information on components traversed
    max_iterations : 1000000
        Flag to stop, in case of any unknown special cases that
        would cause an infinite loop.
    Returns
    ----------
    list of (d,k,comp,rep)
    d : float
        distance to component on B.
    k : float
        k-value of point along segment a.        
    comp : (bool, int)
        component on B 
    rep : (bool, float, float)
        distance representation of component
    """
    # We're going to "walk" along segment a from vertex a to vertex a+1,
    # keeping track of the nearest component on B as we go

    # make sure near component to k=0 is in list
    if not startComp in comp_list:
        comp_list.append(startComp)

    # pre-calculate distance representations and effective intervals
    drs = [hu.distanceRepresentation(A, B, a, c) for c in comp_list]
    eis = [hu.effectiveInterval(A, B, a, c) for c in comp_list]
   
    # intitialize to start of segment a (k=0)
    id = comp_list.index(startComp)
    ei = eis[id]
    k = 0

    # list of (k, d, component, distance representation) along traversal of a
    startRep = hu.distanceRepresentation(A,B,a,startComp)
    nearest_comps = [(startDist,0,startComp,startRep)]

    # in some cases, floating point precision errors will cause 
    # nearest component to have no effective interval.
    # for such cases, set the initial effective interval to [0,1]
    if ei==(float('-inf'),float('-inf')):
        ei = (0,1)

    if verbose:
        print("candidate components: {}".format(len(comp_list)))
        print("\nMoving from {} (ei: {:.3f} - {:.3f}".format(hu.component_label(startComp),ei[0],ei[1]))
    
    # keep going until end of line segment is reached
    numIts = 0
    while id != -1 and numIts < max_iterations:
        numIts += 1
        # move to next component
        if verbose:
            print("initial effective interval: ({:.3f},{:.3f})".format(ei[0],ei[1]))
        new_id,ei,d,k = _updateComponent(A,B,a,comp_list,eis,drs,id,ei,verbose)
        # record components, distance at switch point
        if new_id != -1: # and comp_list[id] != nearest_comps[-1][2]:
            nearest_comps.append((d,k,comp_list[new_id],drs[new_id]))
            # *** CANNOT REMOVE LINE COMPONENTS FROM LIST
            # # remove previous component from list if it is a vertex
            if comp_list[id][0] == False: # component is a vertex
                # delete component from all lists    
                del comp_list[id]
                del drs[id]
                del eis[id]
                # move to next component
                if new_id > id: # we deleted id so subtract one from new_id
                    new_id = new_id - 1
        id = new_id
        
        # report component we're moving to
        if verbose and id != -1:
            comp_label = hu.component_label(comp_list[id])
            print("moving to {} k={:.3f} ei=({:.3f},{:.3f}) d={:.3f}".format(comp_label,k,ei[0],ei[1],d))
    # add in final component at end of segment
    endRep = hu.distanceRepresentation(A, B, a, endComp)
    nearest_comps.append((endDist,1,endComp,endRep))
    if verbose:
        print("moving to {} k={:.3f} d={:.3f}".format(hu.component_label(endComp),1,endDist))

    return nearest_comps

def _updateComponent(A,B,a,comps, eis, drs, prev_id, prev_ei, verbose=False, tol = 0.000001):
    """
    Finds the next component of B while walking along segment a of A. 
    Simultaneously updates the current effective interval

    Parameters
    ----------
    A :     [(x,y),...] list of tuples 
            The coordinates of the main polyline.
    B :     [(x,y),...] list of tuples 
            The coordinates of the other polyline.
    a :     int
            The index of a segment on polyline A.
    comps : [component]
            List of candidate components on B that could be target of Hausdorff distance.
    eis :    list of (float,float)
            K-values of effective intervals of above comps on segment a of A.
    drs :   list of (bool,float,float)
            Distance representations of above comps wrt seg a of A
    prev_id : int
            ID in comps of previous component of polyline B.
    prev_ei : (float,float)
            Previous effective interval on seg a before update
    tol:    float
            Acceptable error tolerance. This is used to avoid floating-point
            precision errors. Adjust at your own risk.

    Returns
    ----------
    (int, [float,float], float)
        compid :        The id wrt comps of the next component of polyline B, 
                        or -1 if no next component is found
        [float,float] :  The updated effective interval of component
        float :          The distance from the switch point between next
                         and previous component to segment a
        float :          The k-value of the switch point
    """
    # NOTES:
    # A key omission in Hangouet seems to be how to handle 3-way intersections
    # in the graph of distance along a segment. These are not so rare at the 
    # start of a segment
    # when the two polylines share some vertices, and in theory they 
    # could happen anywhere.
    
    # To handled botuh three way intersection and floating-point precision errors:
    # (a) define a small "tolerance" (tol)
    # (b) Before switching to a new component at position k, 
    #       > identify all alternate components with switch point less than k + tol
    #       > choose component with shortest distance at k+tol (if in e.i.)
    #       > include current component (i.e. no "switch") in above comparison
    #       > make sure that effective interval starts at >= k + tol


    # *********************
    # NEW METHOD    

    # FIND NEXT SWITCH POINT
    # initialize variables
    cand_k = [2] * len(comps) # list of switch points for each component
    new_comp_id = -1
    new_ei = []
    k_min = 2
    # check all components, looking for one with earliest switch point with previous component
    if verbose:
        comp_label = hu.component_label(comps[prev_id])
        print("   FUNCTION: polyline_hausdorff.update_component")
        print("   prev comp: {}   ei: ({:.3f},{:.3f})".format(comp_label,prev_ei[0],prev_ei[1]))
        print(f"    candidates: {[hu.component_label(c) for c in comps]}")
        dummy = 0
    for cand_id in range(len(comps)):        
        if cand_id != prev_id: # can't return same component!
            if verbose:
                cand_label = hu.component_label(comps[cand_id])
                prev_label = hu.component_label(comps[prev_id])
                print("   testing transition from {} to {}".format(prev_label,cand_label))
            cand_dr = drs[cand_id] #  distance representation of new candidate
            cand_ei = eis[cand_id] #  effective interval of new candidate
            if verbose:
                print('      ei: ({:.3f},{:.3f})'.format(cand_ei[0],cand_ei[1]))
            if len(cand_ei) > 0 and cand_ei != (float('-inf'),float('-inf')): # check if component has an effective interval
                ks = hu.switchPoint(drs[prev_id], cand_dr) # k-values of switch points w/ current component and new candidate
                if verbose:
                    # show calculation of distance representations
                    dr_previous = hu.distanceRepresentation(A, B, a, comps[prev_id])
                    dr_candidate = hu.distanceRepresentation(A, B, a, comps[cand_id])
                    if len(ks) == 0:
                        print("      no crossings found")
                for k in ks:
                    if verbose:
                        if (k < 0 - tol) or (k > 1):
                            print(f"      crossing k ({k}) is not between 0-1")
                        else:
                            print(f"crossing k ({k}) seems okay...")
                    if 0 - tol <= k <= 1: # crossing k is between 0-1
                        if verbose:
                            if (k < prev_ei[0]-tol) or (k > prev_ei[1]+tol):
                                print("      cross point in effective interval of previous component")
                        if prev_ei[0]-tol <= k <= prev_ei[1]+tol:  
                            # Cross point in effective interval of previous component
                            if verbose:
                                if (k < cand_ei[0]-tol) or (k > cand_ei[1]+tol):
                                    print("      Cross point in effective interval of this candidate")
                            if cand_ei[0]-tol <= k <= cand_ei[1] + tol:   
                                # cross point in effective interval of new candidate
                                # note we want to go forward not backward
                                if verbose:
                                    if cand_ei[1] <= prev_ei[0] - tol:
                                        print("      moving backwards")
                                if cand_ei[1] > prev_ei[0] - tol: 
                                    # moving forward    
                                    if k < cand_k[cand_id]:
                                        cand_k[cand_id] = k
                                    if verbose:
                                        if k >= k_min:
                                            print("      k >= k_min")
                                    if k < k_min: # probably this is redundant!
                                        # we might have a new winner
                                        # but need to check distance at tolerance
                                        passes_tolerance_test = True
                                        if k < prev_ei[0] + tol:
                                            # only need to perform test if we have moved less than tolerance
                                            k_check = prev_ei[0] + tol
                                            prev_dist = hu.componentDistance(drs[prev_id], k_check, 1)
                                            new_dist = hu.componentDistance(drs[cand_id], k_check, 1)
                                            if new_dist > prev_dist:
                                                if verbose:
                                                    print("      fails tolerance test!")
                                                passes_tolerance_test = False
                                        if passes_tolerance_test:
                                            if verbose:
                                                print("      PASSED ALL TESTS!!!")
                                            k_min = k
                                            new_comp_id = cand_id
    
    # update effective interval, incrementing by at least tolerance
    comp_label = hu.component_label(comps[new_comp_id])
#     print("   SWITCH POINT: {:.3f} ({})".format(k_min,comp_label))
    cand_ei = eis[new_comp_id]
    new_ei = (max(k_min,prev_ei[0] + tol),cand_ei[1])
    if verbose:
        if new_comp_id != -1:
            print("   winner: {}".format(hu.component_label(comps[new_comp_id])))
        else:
            print("   no switch point found, moving to next component...")

    # If no switch point is found and we haven't reached the end of the segment, 
    # move to adjacent node or vertex
    # *** Hangout says to move to next sequential vertex, 
    #     but actually we have to check both sides
    if new_comp_id == -1 and prev_ei[1] < 1:
        # let's see if this is ever invoked
        comp_label = hu.component_label(comps[prev_id])
        # (print)("   adjacent component search invoked from {}".format(comp_label))
        # get adjacent components
        b = comps[prev_id][1]        
        if comps[prev_id][0] == True: # segment
            adj_comps = [(False,b+i) for i in [0,1] if hu.compValid(len(B),(False,b+i))]
        else: # vertex
            adj_comps = [(True,b+i) for i in [-1,0] if hu.compValid(len(B),(True,b+i))]
        adj_ids = [comps.index(c) for c in adj_comps if c in comps]
        # choose component with lowest effective interval above k    
        if verbose:
            for i in adj_ids:
                comp = comps[i]
                ei = eis[i]
                comp_label = hu.component_label(comp)
                print('   comp: {} ei: ({:.3f},{:.3f})'.format(comp_label,ei[0],ei[1]))
        adj_ids = [i for i in adj_ids if max(eis[i]) > prev_ei[1]] # effective interval above k
        if len(adj_ids) == 0:
            if verbose:
                print("   NO ADJACENT COMPONENTS")
        else:    
            if verbose:
                print("   CHOOSING FROM TWO ADJACENT COMPONENTS")
                print("   prev_ei: ({:.3f},{:.3f})".format(prev_ei[0],prev_ei[1]))
                for i in adj_ids:
                    dr = drs[i] #  distance representation of new candidate
                    ei = eis[i] #  effective interval of new candidate
                    ks = hu.switchPoint(drs[prev_id], dr) # k-values of switch points between current component and new candidate
                    #print('   comp: {} ei: ({:.3f},{:.3f})  ks: {:.3f}'.format(comps[i],ei[0],ei[1],ks))
            
            new_comp_id = min(adj_ids, key = lambda id: min(eis[id])) # lowest
            new_ei = hu.withinUnitInterval(eis[new_comp_id][0],eis[new_comp_id][1])
            # use new component only if we are moving forward along a
            if max(new_ei) < min(prev_ei):
                return (-1,None,None,None)
            
    # check that a new component was found
    if new_comp_id == -1:
        return (-1,None,None,None)
    else:
        # calculate distance at switch point        
        k = new_ei[0]        
        len_a = g.distance(A[a],A[a+1])
        d = hu.componentDistance(drs[new_comp_id], k, len_a)
        if verbose:
            print("   updating to {}".format(hu.component_label(comps[new_comp_id])))
            print("   distance: {:.5f}".format(d))
        return (new_comp_id,new_ei,d,k)
