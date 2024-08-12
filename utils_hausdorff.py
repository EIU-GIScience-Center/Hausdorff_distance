# -*- coding: utf-8 -*-
"""
Hausdorff Utilities (utils_hausdorff.py)
Geometry, vertex and segment components, distance representations,
effective intervals, switch points and other functions required to
calculate the house doff distance. Most of these are taken from
Hangouet 1995, with modifications.

Created on Wed May 26 13:29:02 2021
@author: bjkronenfeld
"""
import math as m
import utils_geom as g
import numpy as np
from rtree import index

def distAcross(comp,dr,interval,len_a,n = 50):
    if comp[0]:
        # component is a segment
        return segDistAcross(dr, interval, len_a)   
    else:
        return vertDistAcross(dr, interval, len_a)

def vertDistAcross(vdr,interval, len_a,n=50):
    """
    returns a list of points defining the distance function from the vertex
    
    Parameters
    ----------
    vdr : vertex distance representation
        As defined based on Hangouet
    interval : (float, float)
        Interval over which to get distances
    len_a : float
        Length of the line segment on which the vertex distance representation 
        is defined
    n : int
        The number of equally spaced distance values to return.
        
    Returns
    ----------
    [(k,d)]
        List of tuples of (k-value along the line segment, distance)
    
    """
    # returns a list of points defining the distance function from the vertex
    r = []
    interval_length = min(interval[1],1)-max(interval[0],0)
    for k in np.arange(max(interval[0],0),min(interval[1],1),interval_length/n):
        d = vertDistance(vdr, k, len_a)
        r.append((k,d))
    return r


def segDistAcross(sdr,interval,len_a,n=50):
    """
    returns a list of points defining the distance function from the segment
    
    Parameters
    ----------
    sdr : segment distance representation
        As defined based on Hangouet
    interval : (float, float)
        Interval over which to get distances
    len_a : float
        Length of the line segment on which the segment distance representation 
        is defined
    n : int
        The number of equally spaced distance values to return.
        
    Returns
    ----------
    [(k,d)]
        List of tuples of (k-value along the line segment, distance)
    
    """

    # returns a list of points defining the distance function from the segment
    r = []
    interval_length = min(interval[1],1)-max(interval[0],0)
    for k in np.arange(max(interval[0],0),min(interval[1],1),interval_length/n):
        d = segDistance(sdr, k, len_a)
        r.append((k,d))
    return r

def componentDistance (dr,k,len_a):
    """
    Calculates the distance from a location on segment a on polyline A to 
    a component of polyline B.
    
    Parameters
    ----------
    A : list of (x,y) tuples
        First polyline.
    dr : distance representation
        distance representation of component of polyline B.
    k : float
        k-value of location along segment a.
    a : int
        ID of segment on polyline A.

    Returns
    -------
    Float
        Distance from location k on segment a to dr component.
    
    # LJ and MT
    """
    
    if dr[0] == True:
        answer = segDistance(dr,k,len_a)
    else:    
        answer = vertDistance(dr,k,len_a)
    return answer

def component_label(comp,verbose = False):
    """for debugging"""
    if comp == None:
        return "None"
    else:
        if verbose:
            if comp[0] == True:
                return f"seg{comp[1]}-{comp[1] + 1}"
            else:
                return f"vrt{comp[1]}"
        else:
            if comp[0] == True:
                return"s{}".format(comp[1])
            else:
                return "v{}".format(comp[1])

def compValid(n_vert,comp):
    """
    Determines if the input component is valid on the input polyline.
    Args:
        polyline : list of tuple
            The input polyline.
        comp : (bool, float)
            The input component
        
    Returns:
        Bool
    
    # BK
    """
    if comp[1] < 0:
        return False
    if comp[0]: # component is a segment
        if comp[1] > n_vert-2:
            return False
    else: # component is a vertex
        if comp[1] > n_vert-1:
            return False
    return True

def distanceRepresentation (A,B,a,bcomp):
    """
    Computes the distance representation of a component of B with respect to segment a 
    Args:
        A,B - the two polylines
        a (int) - segment of A
        bcomp  - componet of B 
        
    Returns:
        distance representation
    
    # LJ and TJ
    """
    # if bcomp is a vertex
    if bcomp[0] == False:
        answer = vertDistRep(A, B, a, bcomp[1])
        return answer
    #if bcomp is a segment
    else:
        answer = segDistRep(A, B, a, bcomp[1])
        return answer


def effectiveInterval(A,B,a,bcomp):
    """
    Determines the range along segment a that is closer to the component bcomp 
    than to either adjacent component.
    Parameters
    ----------
    A : [(x,y),(x,y)...]
        The main polyline.
    B : [(x,y),(x,y)...]
        The other polyline
    a : Integer
        Segment on A
    bcomp : Component
        component of B

    Returns
    -------
    tuple of two k-values

    @author: Folaitan & TJones
    """
    # if statement to determine if bcomp is a segment or vertex
    if bcomp[0] == True:
        sint = segEffectiveInterval(A, B, a, bcomp[1])
        return sint
    elif bcomp[0] == False:
        vint = vertEffectiveInterval(A, B, a, bcomp[1])
        return vint

def segDistance(dr,k,len_a):
    """
    Computes the distance from dr to k on a segment a.
    
    Parameters
    ----------
    dr : (is_seg,sk or None,sin_theta or q-value)
        dist rep of segament B.
        
    k: k-value along segment a.
    
    len_a: length of segment a. 
    Returns
    ----------
    float
        The distance from dr to k on segment a.
        
    @author: megshithakur and Farouk
    """
    # check if segments are parallel
    if dr[1] == None:
        return abs(dr[2] * len_a)
    else: # normal case
        # segDistance formula
        segDist=abs((k-dr[1])*len_a*dr[2])
        return segDist

def segDistRep (A,B,a,b):
    """
    Constructs the component distance representation for segment b with respect to segment a. 
    
    Parameters
    ----------
        A : list of (x,y) tuples
            The first polyline.
        B : list of (x,y) tuples
            The second polyline.
        a : int
            Index of segment on A
        b : int
            Index of segment on B
    
    Returns
    ----------
    tuple (distance representation):
        isSeg : boolean
            True (always True by definition)
        k : float
            k-value of intersection between two segments, or None if segments
            are parallel
        sin_theta : float
            sine of angle between lines through segments, 
            but negative if tangent is negative
            or q-distance between lines if they are parallel.
        k_max : float
            k-value of right hand vertex
    
    # LJ and TJ; k_max added by BK
    """

    qb1 = g.distance_to_line(B[b],A[a],A[a+1],include_sign = True)
    qb2 = g.distance_to_line(B[b+1],A[a],A[a+1],include_sign = True)
    prjb1 = g.project_pt_to_line(B[b],A[a],A[a+1])
    kb1 = g.kvalue(prjb1,A[a],A[a+1])
    prjb2 = g.project_pt_to_line(B[b+1],A[a],A[a+1])
    kb2 = g.kvalue(prjb2,A[a],A[a+1])

    # initiate list for results
    dist_rep = []  
    # Create a boolean statement saying that this is indeed a segment
    dist_rep.append(True)
    # run intersection tool to find the point of intersection between two infinite lines
    x = g.intersection(A[a], A[a+1], B[b], B[b+1])
    # if we've found an intersection, determine sin theta
    if x[0] != None:
        # If there is an intersection, where is it on the x-axis, and then append it
        k_int = g.kvalue(x, A[a], A[a+1]) # convert to a k-value
        dist_rep.append(k_int)
        # Find the angle created by the two lines
        s_rad = g.angle([A[a],A[a+1]],[B[b],B[b+1]])
        # Calculate the sine of the angle
        s = abs(m.sin(s_rad))
        # determine lean of b with respect to a
        # use endpoint furthest from a to determine

        if abs(qb1) > abs(qb2):
            qb, prjb, kb = qb1, prjb1, kb1
        else:
            qb, prjb, kb = qb2, prjb2, kb2
        if (qb > 0) != (kb > k_int):
            # b leans backwards, so we will make sine theta negative
            s = -s
        #Append the sine of the angle
        dist_rep.append(s)
    # Handle case of no intersection
    # If sin theta is zero, there really isn't an intersection (but this 
    # might not be caught above due to floating point precision errors)
    if x[0] == None or dist_rep[2] == 0:
        # reset
        dist_rep=[]
        dist_rep.append(True)
        # Append a no result (no intersection) into the final list
        dist_rep.append(None)
        # Find the distance between the two lines
        q = g.distance_to_line(B[b],A[a],A[a+1])
        # normalize by length of a
        q=q/g.distance(A[a], A[a+1])
        # make negative if b is to the right of the line
        if g.area([A[a],A[a+1],B[b]]) > 0:
            q = -q
        # Append the distance
        dist_rep.append(q)
    # add k-value of right hand of vertex
    
    dist_rep.append(max(kb1,kb2))
    # convert to tuple before returning
    dist_rept = tuple(dist_rep)
    return dist_rept


def withinUnitInterval(a,b):
    """
    Computes the portion of the input interval that is within the interval [0,1]

    Parameters
    ----------
    interval : (float,float)
        An effective interval, not necessarily in sequence.

    Returns
    -------
    The portion of the input interval within [0,1], in sequence from low to high,
    or (-inf,-inf) if the input interval does not overlap the unit interval

    """
    kmin = max(min(a,b),0)
    kmax = min(max(a,b),1)
    if kmin==1 or kmax==0:
        return (float('-inf'),float('-inf'))
    
    else:
        return (kmin,kmax)


def segEffectiveInterval(A,B,a,b, tolerance=0.000001):
    """
    Computes the effective interval of segment b on segment a, that is the interval
    on a for which the interior of segment b is closer than either endpoint
    
    Parameters
    ----------
        A : list of (x,y) tuples
            The first polyline.
        B : list of (x,y) tuples
            The second polyline.
        a : int
            Index of segment on A.
        b : int
            Index of segment on B.
        tolerance : float
            If seg b's k-values on seg a are within this tolerance, 
            the segments will be treated as perpendicular
    Returns:
        tuple
        k1 : float
            First k-value of effective interval. 
        k2 : float
            Second k-value of effective interval.
            
    @author: megshithakur and tannerjones
    """
    a1 = A[a]
    a2 = A[a+1]
    b1 = B[b]
    b2 = B[b+1]
    # check for perpendicular segments
    prjout1 = g.project_pt_to_line(b1, a1, a2)
    prjout2 = g.project_pt_to_line(b2, a1, a2)
    k1 = g.kvalue(prjout1,a1,a2)
    k2 = g.kvalue(prjout2,a1,a2)
    if abs(k1-k2) <= tolerance:
        # segment is perpendicular, so 
        if (g.area([a1,a2,b1]) > 0) == (g.area([a1,a2,b2]) > 0):
            # IF THEY DON'T CROSS effective interval is nothing
            return [float('-inf'),float('-inf')]
        else:
            # IF THEY CROSS effective interval is entire segment 
            return [0,1]
    else:
        # project out from each b vertex to segment a 
        k1 = g.project_out(a1,a2,b1,b2)
        k2 = g.project_out(a1,a2,b2,b1)
        # return them in sequence, bound to range [0,1]
        return withinUnitInterval(k1,k2)


def switchPoint (dr1, dr2, tol = 0.00000001, verbose = False):
    """
    Determines point along a where the nearest point on B switches from the 
    component represented by dr1 to the component represented by dr2.

    Parameters
    ----------
    dr1 : dr1
        Distance representation 1
    dr2 : dr2
        Distance representation 2
    b1 : component
        first component of B
    b2 : component
        second component of B
    B : coordinates of polyline B

    Returns
    -------
    [float]	list of (at most 2) k-values

    """
    # input: segment, segment - return Seg Seg switch point
    if dr1[0] == True and dr2[0] == True:
        segseg = segSegSwitchPoint(dr1, dr2)
        return segseg
    # input: vertex, segment - return vert seg switch point
    elif dr1[0] == False and dr2[0] == True:
        vertseg = vertSegSwitchPoint(dr1, dr2, from_vertex = True, tol = tol)
        return vertseg
    # input: segment, vertex - return vert seg switch point
    elif dr1[0] == True and dr2[0] == False:
        vertseg = vertSegSwitchPoint(dr2, dr1, from_vertex = False, tol = tol)
        return vertseg
    # input: vertex, vertex - return vert vert switch point
    elif dr1[0] == False and dr2[0] == False:
        vertvert = vertVertSwitchPoint(dr1, dr2)
        return vertvert


def vertDistance(dr, k, len_a):
    """
    Computes the distance from the vertex represented by dr to the location k 
    on the segment a that the distance representation was constructed from.
    Parameters
    ----------
    dr : (is_seg, vk,q)
        distance representation of a vertex of B, being a tuple of three values
    k : the k-value along segment a
    len_a : the length of segment a

    Returns
    -------
    float
        the distance from point k on segment a to the vertex of B represented 
        by dr
        
    @author: Folaitan & @Ljansen
    """
    d = ((k-dr[1])**2) + (dr[2]**2)
    d = m.sqrt(d)
    d = len_a * d
    return abs(d)


def vertDistRep(A,B,a,b):
    """
    Constructs the Distance Representation for vertex b with respect to segment a.
    
    Parameters
    ----------
        A : list of (x,y) tuples
            The first polyline.
        B : list of (x,y) tuples
            The second polyline.
        a : int
            Index of segment on A
        b : int
            Index of vertex on B
    Returns
    ----------
    tuple (distance representation): 
        is_seg : bool 
            False(always False by definition).
        k : float 
            k-value of the location of the perpendicular projection of b onto the line through a.
        q : float 
            q-distance from b to the line through a, negative if b is right of a
            
    @author: MT and TJ
    """
    # initiate list for results
    FinalList = []
    #Create a boolean statement saying that this is not a segment
    FinalList.append(False)
    #run project to line tool to find perpendicular intersection point between point b and segment a 
    p = g.project_pt_to_line(B[b], A[a], A[a+1])
    k = g.kvalue(p, A[a], A[a+1])
    #append result to results list
    FinalList.append(k)
    #run distance tool to find distance between the new point (k) and b
    q = g.distance(p, B[b])
    q = q/g.distance(A[a],A[a+1])
    # make negative if b is to the right of the line
    if g.area([A[a],A[a+1],B[b]]) > 0:
        q = -q
    #append the result
    FinalList.append(q)
    FinalList = tuple(FinalList) 
    return FinalList


def vertEffectiveInterval(A,B,a,b,tolerance=0.000001):
    """
    Computes the effective interval of vertex b on segment a, that is the interval
    on a for which vertex b is closer than either adjacent segment
    
    Parameters
    ----------
        A : list of (x,y) tuples
            The first polyline.
        B : list of (x,y) tuples
            The second polyline.
        a : int
            Index of segment on A
        b : int
            Index of vertex on B
        tolerance : float
            If b segs' k-values on seg a are within this tolerance, 
            the segments will be treated as perpendicular
    Returns:
    ----------
        tuple of two floats, with (-inf,-inf) representing no effective
        interval.
        
    
    @author: megshithakur and tannerJones
    """
    # get coordinates
    a1, a2 = A[a],A[a+1]

    
    # handle cases where b is an end of the polyline
    if b==0 or b==len(B)-1:
        # get neighbor segment id
        if b==0:
            nb = 1
        else:
            nb = len(B)-2
        # get coordinates of b vertex and neighbor
        bc,bnb = B[b],B[nb] 
        
        # determine projections onto segment a
        prjc = g.project_pt_to_line(bc,a1,a2)
        prjnb = g.project_pt_to_line(bnb,a1,a2)
        # determine k-values
        kc = g.kvalue(prjc, a1, a2)
        knb = g.kvalue(prjnb, a1, a2)

        # check for perpendicularity
        if abs(kc-knb) <= tolerance:
            # calculate distances from each B vertex to segment a
            distc = g.distance(bc,prjc)
            distnb = g.distance(bnb,prjnb)
            # is vertex b on segment a?
            if distc < tolerance:
                return (0,1)
            # is the neighboring vertex on segment a?
            elif distnb < tolerance:
                return (float('-inf'),float('-inf'))
            # are the two vertices on opposite sides of segment a?
            else: 
                # calculate areas
                areac = g.area([a1,a2,bc,a1])
                areanb = g.area([a1,a2,bnb,a1])
                # are the two vertices on opposite sides of segment a?
                if (areac > 0) != (areanb > 0):
                    return (float('-inf'),float('-inf'))
                # is the center vertex closer to segment b than the neighbor
                elif distc < distnb:
                    return [0,1]
                else:
                    return (float('-inf'),float('-inf'))
        else:
            # project points out from b to segment a
            kbout = g.project_out(a1, a2, bc, bnb)
            # return interval from kbout to end of segment a opposite neighbor
            if knb < kc:
                return withinUnitInterval(kbout,1)
            else:
                return withinUnitInterval(0,kbout)
    else: # vertex b has two neighbors     
        # get coordinates of B vertices
        bp,bc,bn = B[b-1],B[b],B[b+1] # prev, current, next vertices on B
        # project each vertex of B onto segment A
        prjp = g.project_pt_to_line(bp,a1,a2)
        prjc = g.project_pt_to_line(bc,a1,a2)
        prjn = g.project_pt_to_line(bn,a1,a2)
        # get k-values of projections of B vertices onto a
        kp = g.kvalue(prjp, a1, a2)
        kc = g.kvalue(prjc, a1, a2)
        kn = g.kvalue(prjn, a1, a2)    
        # check for perpendicular segments
        prev_perp = abs(kc-kp) <= tolerance
        next_perp = abs(kn-kp) <= tolerance
        if prev_perp and next_perp: # both perpendicular
            # get distances and areas to all vertices
            dcur = g.distance(bc,prjc)
            dprev = g.distance(bp,prjp)
            dnext = g.distance(bp,prjn)
            areacur = g.area([a1,a2,bc,a1])
            areaprev = g.area([a1,a2,bp,a1])
            areanext = g.area([a1,a2,bn,a1])
            # vertex has effective interval only if it is closest to segment a
            # and all vertices are on same side of segment a
            if (areaprev > 0) == (areanext > 0) and (areacur > 0) == (areanext > 0) and dcur < dprev and dcur < dnext:
                return (0,-1)
            else:
                return (float('-inf'),float('-inf'))
        elif prev_perp or next_perp: # one perpendicular
            # get coordinates and k-values of vertex on perpendicular segment, other vertex
            if prev_perp:
                bperp = bp
                bother = bn
                kother = kn
                prjperp = prjp
            else:
                bperp = bn
                bother = bp
                kother = kp
                prjperp = prjn            
            # get distances and areas of each vertex on perpendicular segment
            dperp = g.distance(bperp,prjperp)
            dcur = g.distance(bc,prjc)
            areaperp = g.area([a1,a2,bperp,a1])                    
            areacur = g.area([a1,a2,bc,a1])                    
            # Is perpendicular segment on one side of A and current vertex is closer?
            if (areaperp > 0) == (areacur > 0) and dcur < dperp:
                # project out from other segment
                otherprjoutk = g.project_out(a1,a2,bc,bother)
                # interval is from end of line opposite other to k-value of other
                if kother > kc:
                    return (0,otherprjoutk)
                else:
                    return (otherprjoutk,1)
            else:
                return (float('-inf'),float('-inf'))
        else: # neither perpendicular
            # get k-values of projections of vertex out from each B segment onto segment a            
            kcpout = g.project_out(a1, a2, bc, bp) 
            kcnout = g.project_out(a1, a2, bc, bn)
            # check sides of neighboring vertices with respect to b
            if kp < kc and kn < kc: # both neighbors left of b
                maxk = max(kcpout,kcnout)
                return withinUnitInterval(maxk,1)
            elif kp > kc and kn > kc: # both neighbors right of b            
                mink = min(kcpout,kcnout)
                return withinUnitInterval(0,mink)
            else: # nieghbors on either side of b
                # determine min and max k-values of b based on positions of
                # previous and next vertices
                if kp < kc: # previous neighbor left
                    mink = kcpout
                    maxk = kcnout
                else: # previous neighbor right
                    mink = kcnout
                    maxk = kcpout
                if mink <= maxk: 
                    return withinUnitInterval(mink,maxk)
                else:
                    return (float('-inf'),float('-inf'))


def segSegSwitchPoint(dr1, dr2, tol = 0.00000001):
    """
    Determines the k-value where segment A crosses from 
    the Voronoi cell of the first segment
    into the Voronoi cell of the second segment

    Parameters
    ----------
    dr1 : seg_p1
        distance representation of a segment on b 
    dr2 : seg_p2
        distance representation of another segment on b 

    Returns
    -------
    List of 1 float representing k-value of the switch points
    
    @author: Folaitan & Ljansen; directional component added by BK
    """
    k = []
    if dr1[1] == None and dr2[1] == None: # both segments parallel to a
        return []
    elif dr1[1] == None: # first segment parallel to a
        q = dr1[2] # distance from a to parallel segment
        if abs(q) > tol:
            k_out_1 = dr2[1] + (dr1[2]/dr2[2])
            k_out_2 = dr2[1] - (dr1[2]/dr2[2])
            k.append(min(k_out_1,k_out_2))
        else:
            return []
    elif dr2[1] == None: # second segment parallel to a
        k_out_1 = dr1[1] + ((dr2[2])/dr1[2])
        k_out_2 = dr1[1] - ((dr2[2])/dr1[2])
        k.append(max(k_out_1,k_out_2))
    else: # neither segment parallel to a
        # create more readable variables
        k1 = dr1[1]
        k2 = dr2[1]
        s1 = dr1[2]
        s2 = dr2[2]
        # calculate alpha parameter
        a = s2/s1
        if a == 1: # two segments are parallel
            # first solution is not valid
            # return second solution only if k2 > k1
            if k2 > k1:
                k_out_2 = (k1+a*k2)/(1+a)
                k.append(k_out_2)
            else:
                return []
        elif a == -1: # two segments are mirror images of each other
            # second solution is not valid
            # return first solution only if k2 > k1
            if k2 > k1:
                k_out_1 = (k1-a*k2)/(1-a)
                k.append(k_out_1)
            else:
                return []
        else: # there are two solutions
            k_out_1 = (k1-a*k2)/(1-a)
            k_out_2 = (k1+a*k2)/(1+a)
            # select solution based on which statement is further right
            if abs(s2) < abs(s1):
                # b2 is closer to parallel with segment a
                k.append(max(k_out_1,k_out_2))
            else:
                # b1 is closer to parallel the segment a
                k.append(min(k_out_1,k_out_2))
    # return values in ascending order, for consistency
    return sorted(k)


def vertSegSwitchPoint(vdr, sdr, from_vertex = True, tol = 0.00000001):
    """
    Determines the k-value of the location along segment a that is equidistant 
    to the input vertex and input segment, i.e. the location at which the nearest 
    component “switches” from the first component to the second.
    Parameters
    ----------
    vdr : ver_rep : (False, k, q)
        distance representation of a vertex on B
    sdr : seg_rep : (True, k, sin_theta) or (True, none, q)
        distance representation of a segment on B
    ks_max : float
        The maximum K value of the input segment
    from_vertex : Boolean
        true if we are looking to go from the vertex to the segment
        false if we are looking to go from the segment to the vertex

    Returns
    -------
    [float]
    List containing 0 or 1 floats representing the k-values of 
    the switch point if it exists
    
    @author: Folaitan & Ljansen; directionality introduced by BK
    """
    verbose = False
    # get distance representation values into more readable variables
    k_vert = vdr[1]
    q_vert = vdr[2]
    k_seg = sdr[1]
    sin_theta = sdr[2]
    ks_max = sdr[3]
    
    # three cases: segment b perpendicular, parallel, or oblique to A
    if abs(sin_theta) == 1: 
        # b segment is perpendicular to A
        # only return a crossing if second component is right of first        
        if k_seg == k_vert:
            # point and line have same k-value, so no switching
            return []
        elif from_vertex != (k_vert > k_seg):
            # going from vertex to segment but vertex is to the right
            # or else going from segment to vertex but segment is to the right
            # so no switching
            return []
        else:
            # the component we want to go to is to the right
            # so make the switch
            numerator = ((k_seg**2)-(k_vert**2)-(q_vert**2))
            denominator = (2*k_seg) - (2*k_vert)
            vertSegPoint = numerator/denominator
            return [vertSegPoint]
    elif k_seg == None: 
        # b segment is parallel to a
        if verbose:
            print("Segment is parallel to a...!!!!!!!")
        q_seg = sin_theta # distance representation value is q not sin_theta
        a = 1
        b = -2*k_vert
        c = k_vert**2 + q_vert**2 - q_seg**2
    else: # normal case
        if verbose:
            print("Normal case!!!!!!!")
        a = (sin_theta**2)-1
        b = (2*k_vert)-((2*k_seg)*(sin_theta**2))
        c = ((k_seg**2)*(sin_theta**2))-(k_vert**2)-(q_vert**2)
    inside_root = (b**2)-(4*a*c)
    # catch floating point precision issues:
    if -tol < inside_root < tol:
        if verbose:
            print("Setting inside routes to zero due to floating point precision...")
        inside_root = 0                
    # check number of roots
    if inside_root < 0:
        if verbose:
            print("Inside route is less than zero!!!!!!!")
        # if inside_root is less than zero, quadratic formula has no 
        # solution and there is no switch point to return
        return []
    elif inside_root == 0:
        # if it equals zero, there is one switch point
        # we have either a tangency or point is on line
        try:
            point_on_line = True
            if k_seg == None: 
                if verbose:
                    print("no k_seg!!!!!!!")
                # segment is parallel to a
                # just check if segment and vertex are on same side of a
                if (q_vert > 0) != (q_seg > 0):
                    point_on_line = False
            ### NEW METHOD
            elif abs(q_vert - (k_vert - k_seg) * (sin_theta / (1-sin_theta ** 2) ** 0.5)) > 2 * tol:
                point_on_line = False
            ### OLD METHOD
            # elif (k_seg > k_vert) == ((q_vert > 0) == (sin_theta > 0)):
            #     print("really complicated triple equality!!!!!!!")
            #     print(f"k_seg: {k_seg}")
            #     # the above returns True if 0 or 2 of the terms are true
            #     # in which case segment B is tangent to parabolic Voronoi edge
            #     point_on_line = False
            ### END OLD METHOD
                
                
            if point_on_line:
                # vertex is on segment; need to check which is right and
                # make sure we get the vertex
                if verbose:
                    print("point is online!!!!")
                if from_vertex and (ks_max - tol > k_vert):
                    return [-b/(2*a)]
                elif (from_vertex == False) and (k_vert >= ks_max - tol):
                    return [-b/(2*a)]
                else:
                    return []
            else:
                # no need to switch
                if verbose:
                    print("point is not on line!!!!!!!")
                return []

        except Exception as e:
            print(e)
            print("k_seg: {}".format(k_seg))
            print("k_vert: {}".format(k_vert))
            print("q_vert: {}".format(q_vert))
            print("sin_theta: {}".format(sin_theta))
            
    if inside_root > 0:
        # we have two solutions,  answer depends on which way we are going        
        solution_1 = (-b - (m.sqrt(inside_root)))/(2*a)
        solution_2 = (-b + (m.sqrt(inside_root)))/(2*a)
        if from_vertex:
            # pick larger solution
            return [max(solution_1, solution_2)]
        else:
            return [min(solution_1, solution_2)]


            
def vertVertSwitchPoint (dr1,dr2):
    """
    Determines the k-value of the location along segment a that crosses from
    the Voronoi cell around v1 into the Voronoi cell around v2

    Parameters
    ----------
    dr1 : tuple, distance representation of a vertex of B
    dr2 : tuple, distance representation of a different vertex of B

    Returns : [k] - list containing one float representing the k-value of the switch point
    
    # MT and LJ; directionality introduced by BK
    """
    # return crossing only if v2 is to the right of v1 (i.e. has higher k-valuep)
    if dr2[1] > dr1[1]:
        #calculates k value for switch point
        k = [((dr2[2]**2 - dr1[2]**2) + (dr2[1]**2 - dr1[1]**2)) / (2 * (dr2[1]-dr1[1]))] 
    else:
        #k is an empty list
        k = []
    return k

def b_box(B,b):
    """ Creates a box around segment b"""
    left = min(B[b][0],B[b+1][0])
    bottom = min(B[b][1],B[b+1][1])
    right = max(B[b][0],B[b+1][0])
    top = max(B[b][1],B[b+1][1])
    return (left,bottom,right,top)    

def seg_idx(B):
    """
    Creates an r-tree index of the segments in polyline B.

    Parameters
    ----------
    B : [(x,y),...] list of tuples 
        A polyline.

    Returns
    -------
    r : r-tree index
        Index for searching.

    """
    # use stream loading, as in 
    # https://rtree.readthedocs.io/en/latest/performance.html
    boxes = [b_box(B,b) for b in range(len(B)-1)]
    def gen_func():
        for i, box in enumerate(boxes):
            yield (i,box,i)
    r = index.Index(gen_func())
    return r


def candidateComponents(A,B,a,nl1,nl2,d1,d2,B_seg_idx,brute_force = False):
    """
    Identifies all components of B that could be the target of the Hausdorff 
    distance from segment a on A.

    Parameters
    ----------
    A : [(x,y),...] list of tuples 
        The coordinates of the main polyline.
    B : [(x,y),...] list of tuples 
        The coordinates of the other polyline.
    a : int
        The index of a segment on polyline A.
    d1 : float
        Distance from vertex a to nearest component on B.
    d2 : float
        Distance from vertex a+1 to nearest component on B.
    B_seg_idx : rtree index
        Index of segments on B.
        
    Returns
    ----------
    [(bool, int)]
        A list of candidate components on B.
    """
    # brute_force = True

    if brute_force:
        # OLD
        # for now, simply return a list of all components of B
        result = []
        # get vertex components
        for i in range(len(B)):
            result.append((False,i))
        for i in range(len(B)-1):
            result.append((True,i))
        return result

    else:
        # NEW - USE R-TREE INDEX FOR EFFICIENCY
        result = []
        # radii of three circles
        r1 = d1 # circle centered on A[a]
        r2 = d2 # circle centered on A[a+1]
        # Center two near locations on B
        bc = ((nl1[0] + nl2[0])/2, (nl1[1] + nl2[1])/2)
        # Slope of line between two near locations on B
        bdx = nl2[0] - nl1[0]
        bdy = nl2[1] - nl1[1]
        # Two points on perpendicular bisector of above line
        perp_1 = (bc[0] + bdy, bc[1] - bdx)
        perp_2 = (bc[0] - bdy, bc[1] + bdx)
        # Point on segment a equidistant to two near locations
        # calculated as intersection between perpendicular bisector and seg a
        anchor = g.intersection(perp_1, perp_2, A[a], A[a+1])
        if anchor == (None,None):
            # if perp bisector is parallel to seg a, anchor is w
            # on line between near locations
            # first check if near locations are the same
            if nl1 == nl2:
                # project point to linegeometry utilities
                anchor = g.project_pt_to_line(nl1,A[a],A[a+1])
            else:
                anchor = g.intersection(nl1, nl2, A[a], A[a+1])
        # third radius is from anchor

        rc = g.distance(anchor,nl1)
        
        # get bounds of search
        xmin = min(A[a][0] - r1, A[a+1][0] - r2, anchor[0] - rc)
        xmax = max(A[a][0] + r1, A[a+1][0] + r2, anchor[0] + rc)
        ymin = min(A[a][1] - r1, A[a+1][1] - r2, anchor[1] - rc)
        ymax = max(A[a][1] + r1, A[a+1][1] + r2, anchor[1] + rc)

        # new method: use anchor only
        xmin = anchor[0] - rc
        xmax = anchor[0] + rc
        ymin = anchor[1] - rc
        ymax = anchor[1] + rc

        # search for candidate components using index        
        square = (xmin,ymin,xmax,ymax)       
        segs = list(B_seg_idx.intersection(square))
        verts = set()
        [verts.add(s) for s in segs]
        [verts.add(s+1) for s in segs]
        verts = list(verts)
        for s in segs:
            result.append((True,s))
        for v in verts:
            result.append((False,v))
        return result


def nearSegment(A,B,a, B_seg_idx):
    """
    Identifies the segment of of B nearest to vertex a on A.

    Parameters
    ----------
    A : [(x,y),...] list of tuples 
        The coordinates of the main polyline.
    B : [(x,y),...] list of tuples 
        The coordinates of the other polyline.
    a : int
        The index of a vertex on polyline A.
    B_seg_idx : rtree index
        Index of B segments for faster searching.
    Returns
    ----------
    int
        The index of the nearest segment on B.
    """
    # for now, this will be coded with a "brute force" method, checking the 
    # distance from every segment of B to vertex a of A
    # later this should be updated to use an indexing structure such as 
    # an r-tree for computational efficiency
    # initialize to first segment

    brute_force = False

    if brute_force:
        # OLD
        min_index = 0
        min_d = g.distance_to_segment(A[a], B[0], B[1])
        # check other segments
        for b in range(1,len(B)-1):
            d = g.distance_to_segment(A[a],B[b],B[b+1])
            if d < min_d:
                min_d = d
                min_index = b
        return min_index

    else:    
        # NEW - USE R-TREE FOR EFFICIENCY
        min_index = 0
        abox = (A[a][0],A[a][1],A[a][0],A[a][1])
        # get index of "first guess" segment with nearest bounding box
        fg = list(B_seg_idx.nearest(abox,1))[0]
    
        # get distance fgd to first guess
        fgd = g.distance_to_segment(A[a],B[fg],B[fg+1])
    
        # get all segments in square of radius fgd around a
        square = (A[a][0]-fgd,A[a][1]-fgd,A[a][0]+fgd,A[a][1]+fgd)
        cands = list(B_seg_idx.intersection(square))
        
        # loop through segments to find closest
        min_d = fgd
        min_index = fg
        # check other segments
        for b in cands:
            d = g.distance_to_segment(A[a],B[b],B[b+1])
            if d < min_d:
                min_d = d
                min_index = b    
        return min_index

    
def nearComponent(A,B,a,b):
    """
    Among segment b, vertex b and vertex b+1, determines which component is nearest to vertex a.

    Parameters
    ----------
    A : [(x,y),...] list of tuples 
        The coordinates of the main polyline.
    B : [(x,y),...] list of tuples 
        The coordinates of the other polyline.
    a : int
        The index of a vertex on polyline A.
    b : int
        The index of a segment on polyline B.
    Returns
    ----------
    (component, (float,float), float)
        The nearest component, pt and distance from vertex a of B
    """
    # # get distances from each component
    # d_seg = g.distance_to_segment(A[a],B[b],B[b+1])
    # d_vert1 = g.distance(A[a],B[b])
    # d_vert2 = g.distance(A[a], B[b+1])
    # # if distance to d_seg is strictly the minimum, return the segment
    # if d_seg == min(d_seg,d_vert1,d_vert2):
    #     return ((True,b),d_seg)
    # elif d_vert1 < d_vert2:
    #     return((False,b),d_vert1)
    # else:
    #     return ((False,b+1),d_vert2)
    
    # calculate k-value of projection of vertex a onto segment b
    prj = g.project_pt_to_line(A[a], B[b], B[b+1])
    k = g.kvalue(prj, B[b], B[b+1])
    if k <= 0: # return vertex b
        d = g.distance(A[a],B[b])
        return ((False,b),B[b],d)
    elif k >= 1: # return vertex b+1
        d = g.distance(A[a],B[b+1])
        return ((False,b+1),B[b+1],d)
    else: # return segment b
        d = g.distance(A[a],prj)
        return ((True,b),prj,d)        


def nearLoc(srcloc,trgline,trgcomp):
    """
    Finds the location on the target component nearest to the source
    location.

    Parameters
    ----------
    srcloc : (float,float)
        The coordinates of the source location.
    trgline : [(x,y),...] list of tuples 
        The coordinates of the target polyline.
    trgcomp : (bool,int)
        The target component.
    
    Returns
    ----------
    (float,float)
        The nearest location on the target component to the source location.
    """
    if trgcomp[0] == False: # target is a vertex
        return trgline[trgcomp[1]]
    else: # target is a segment
        trgstart = trgline[trgcomp[1]]
        trgend = trgline[trgcomp[1]+1]
        srcprj = g.project_pt_to_line(srcloc, trgstart,trgend)
        k = g.kvalue(srcprj,trgstart,trgend)
        if k <= 0:
            return trgstart
        elif k >=1:
            return trgend
        else:
            x = trgstart[0] + k * (trgend[0]-trgstart[0])
            y = trgstart[1] + k * (trgend[1]-trgstart[1])
            return (x,y)
        
def checkSegment(c1,c2):
    """
    Determines whether or not it is necessary to further process a segment
    given that its endpoints have been processed already.

    Parameters
    ----------
    c1 : component
        The component of B closest to the first vertex of a segment of A.
    c2 : component
        The component of B closest to the second vertex of a segment of A.
    Returns
    ----------
    (component, float)
        The nearest component of B along with its distance from vertex a.
    """
    # No need to check a segment if either:
    # c1 and c2 are the same component, or
    # c1 and c2 are consective vertices
    if c1==c2:
        return False
    elif c1[0] == False and c2[0] == False and abs(c1[1]-c2[1])==1:
        return False
    else:
        return True
