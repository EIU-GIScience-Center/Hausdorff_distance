"""
-------------------------------------------------------------------------------
 Name:        geom utils
 Purpose:     various geometry functions
 License:     MIT License
 Updated:     April 2021
 Notes:       Docstrings follow numpy format 
              https://numpydoc.readthedocs.io/en/latest/format.html#docstring-standard
-------------------------------------------------------------------------------
"""

import math as __m
from decimal import Decimal



def angle(v1,v2):
    """
    Computes the angle from vector v1 to vector v2, measured in radians. 

    Parameters
    ----------
    v1 : [(x,y),(x,y)] list of two tuples 
        A vector between two points.
    v2 : [(x,y),(x,y] list of two tuples 
        A vector between two points. Ordinarily the first points of v1 and v2
        will be the same, but this is not necessary.
    Returns
    ----------
    float
        angle between the two segments defined by the input
    """
    # code by BJK following method in 
    #   https://scicomp.stackexchange.com/questions/27689/numerically-stable-way-of-computing-angles-between-vectors
    # which is measured to be numerically stable to approx. 15 decimal places
    
    # move vectors to origin
    norm1=(v1[1][0]-v1[0][0],v1[1][1]-v1[0][1])
    norm2=(v2[1][0]-v2[0][0],v2[1][1]-v2[0][1])
    # get distances of triangle formed
    a=distance((0,0),norm1)
    b=distance((0,0),norm2)
    c=distance(norm1,norm2)
    #a=distance(v1[0],v1[1])
    #b=distance(v2[0],v2[1])
    #c=distance(a,b)

    # calculate mu
    if b<0 or c<0: # this should never happen
        return None
    if b >= c:
        mu = c-(a-b)
    else:
        mu = b-(a-c)
    # calculate numerator and denominator of ratio in formula
    numerator=((a-b)+c)*mu
    denominator=(a+(b+c))*((a-c)+b)
    
    if denominator == 0:
        # if denominator is zero, that means a and c are same length and b is zero
        # so angle is 180 degrees
        return __m.pi
    else:
        # if numerator is less than zero, it is due to floating point precision
        # it should be very close to zero, and the angle is zero
        if -0.000001 < numerator/denominator < 0:
            return 0
        else:
            # otherwise, calculate formula
            half_tangent = __m.sqrt(numerator/denominator)
            theta = 2 * __m.atan(half_tangent)    
            theta = abs(theta)
            if theta > __m.pi:
                theta = 2*__m.pi-theta
            return theta 


def area(pts,absolute=False):
    """    
    Computes the clockwise area of the polygon defined by the points.

    Parameters
    ----------
    pts : list of (x,y) tuples 
        The coordinates of a polygon. 
    absolute : bool
        Flag indicating if return value must be a positive number.
    Returns
    ----------
    float
        The area of the input polygon.
    """
    # uses shoelace formula (surveyor's formula)
    # check to see if polygon closes; if not, then force it to close
    if pts[len(pts)-1] != pts[0]:
        pts.append(pts[0])
    # calculate areas of each trapezoid
    a=[(pts[i+1][0]-pts[i][0])*(pts[i][1]+pts[i+1][1]) for i in range(len(pts)-1)]
    # sum areas and divide by 2
    A=sum(a)/2
    # handle request for absolute value
    if absolute:
        return abs(A)
    else:
        return A


def densify(poly,spacing):
    """
    Inserts vertices into each segment at maximum spacing distance apart

    Parameters
    ----------
    poly : list of (x,y) tuples
        Coordinates of polyline.
    spacing : float
        Maximum spacing between vertices.

    Returns
    -------
    Densified polyline.

    """
    r = [poly[0]]
    for v in poly[1:]:
        prev_v = r[-1]
        d = distance(prev_v,v)
        nv = __m.ceil(d/spacing) - 1
        dx = (v[0] - prev_v[0]) / (nv + 1)
        dy = (v[1] - prev_v[1]) / (nv + 1)
        for i in range(nv):
            x = prev_v[0] + dx * (i + 1)
            y = prev_v[1] + dy * (i + 1)
            r.append((x,y))
        r.append(v)
    return r

#
def distance(a,b):
    """ 
    Computes the distance between two points.
    
    Parameters
    ----------
    a : (x,y)
        Coordinates of a point.
    b : (x,y)
        Coordinates of another point.
    Returns
    ----------
    float
        The Euclidean distance between the two input points.
        
    """
    return  ((b[0]-a[0])**2 + (b[1]-a[1])**2)**0.5
    # return  __m.sqrt((b[0]-a[0])**2 + (b[1]-a[1])**2)


def project_pt_to_line(p,a,b):
    """
    Finds the point on the line through a and b closest to p.
    
    Parameters
    ----------
    p : (x,y)
        Coordinates of a point.
    a : (x,y)
        Coordinates of a point on a line.
    b : (x,y)
        Coordinates of another point on a line.
    Returns
    ----------
    (x,y)
        The point on the infinite line through a & b that is closest to p.
    """
    # code by BJK
    
    # if line is degenerate, this can't be computed
    if a==b:
        return None
    # get line length
    lineD = distance(a,b)
    # get area of triangle formed by three points
    A = area([a,b,p,a])
    # vector distance, which can be negative
    vecD = 2*A/lineD
    # line vector differentials
    lineVecX = b[0]-a[0]
    lineVecY = b[1]-a[1]
    # vector to result point is clockwise perpendicular to line
    rX = -lineVecY
    rY = lineVecX
    # scale by ratio of vector length to baseline length
    scalar = vecD/lineD
    rX = rX*scalar
    rY = rY*scalar
    # move from input point along vector    
    outX = p[0] + rX
    outY = p[1] + rY
    # return as an (x,y) tuple
    return (outX,outY)


def distance_to_line(p,a,b, include_sign = False):
    """
    Computes the perpendicular distance from a point to an infinite line.
    
    Parameters
    ----------
    p : (x,y)
        Coordinates of a point.
    a : (x,y)
        Coordinates of a point on a line.
    b : (x,y)
        Coordinates of another point on a line.
    include_sign : Boolean
        if True, return value will be positive if point is to the left of
        the line, negative if it is to the right
    Returns
    ----------
    float
        The Euclidean distance from p to the infinite line through a & b.
    """
    # code by BJK
    # area of triangle formed between point and line segment
    trianglearea = -area([a,b,p])
    if include_sign == False:
        trianglearea = abs(trianglearea)
    # length of line segment
    line_length=distance(a,b)
    # make sure line segment has a length
    if line_length==0: 
        # a & b are the same, so just calculate distance between points
        return distance(p,a) 
    else: 
        # the distance we want is the height of the triangle
        # area is 1/2 base x height so height is 2*area/base
        return 2*trianglearea/line_length


def is_monotonic(a,b,c):
    """
    Determines if the three input numbers are in sequence, either low-to-high
    or high-to-low (with ties acceptable).
    
    Parameters
    ----------
    a : float
    b : float
    c : float
        
    Returns
    ----------
    boolean
        
    """
    # if b is equal to either a or c then the sequence is monotonic
    if b == a:
        return True
    elif b == c:
        return True
    # otherwise, a and c must be on different sides of b
    elif (a > b) != (c > b):
        return True
    else:
        return False    

    
def distance_to_segment(p,a,b):
    """
    Computes the perpendicular distance from a point to a finite line segment.
    
    Parameters
    ----------
    p : (x,y)
        Coordinates of a point.
    a : (x,y)
        Coordinates of a point on a line.
    b : (x,y)
        Coordinates of another point on a line.
    Returns
    ----------
    float
        The Euclidean distance from p to the finite line segment through a & b.
    """
    # project point onto segment
    p2 = project_pt_to_line(p,a,b)
    # determine if projection point is on segment or not
    if is_monotonic(a[0],p2[0],b[0]) and is_monotonic(a[1],p2[1],b[1]):
        # if so, use distance from p to its projection
        return distance(p,p2)
    else:
        # otherwise, use distance to nearest endpoint 
        return min(distance(p,a),distance(p,b))
    
    
def intersection(A,B,C,D,infinite=True,ultra_precise=False):
    """ 
    Computes the point of intersection between two lines.

    Parameters
    ----------
    A : (x,y)
        Coordinates of a point on the first line.
    B : (x,y)
        Coordinates of another point on the first line.
    C : (x,y)
        Coordinates of a point on the second line.
    D : (x,y)
        Coordinates of another point on the second line.
    infinite : bool
        If false, will only look for an intersection on the line segments.
    ultra_precise : bool
        If true, will use decimals instead of floats (much slower but more precise)
    Returns
    ----------
    (x,y)
        The coordinates of the intersection point between AB and CD, or 
        (None,None) if no intersection point is found
    """
    # code by BJK, following Stephen Wise (2013), GIS Fundamentals, pp. 48-9
    if ultra_precise: 
        # use decimal representations (much slower but precise beyond 8 decimal places)
        A=(Decimal(A[0]),Decimal(A[1])) 
        B=(Decimal(B[0]),Decimal(B[1])) 
        C=(Decimal(C[0]),Decimal(C[1])) 
        D=(Decimal(D[0]),Decimal(D[1])) 
    if A[0]==B[0]:
        if C[0]==D[0]:
            xp,yp=None,None # lines are parallel
        else: # first line vertical
            b2=(D[1]-C[1])/(D[0]-C[0])
            a2=C[1]-b2*C[0]
            xp = A[0]
            yp = a2+b2*xp
    else:
        if C[0]==D[0]: # second line vertical
            b1=(B[1]-A[1])/(B[0]-A[0])
            a1=A[1]-b1*A[0]
            xp = C[0]
            yp = a1+b1*xp
        else: # neither line vertical
            b1=(B[1]-A[1])/(B[0]-A[0])
            b2=(D[1]-C[1])/(D[0]-C[0])
            a1=A[1]-b1*A[0]
            a2=C[1]-b2*C[0]
            if b1==b2:
                xp,yp = None,None # lines are parallel
            else:
                xp = -(a1-a2)/(b1-b2)
                yp = a1+b1*xp
    # test whether intersection point falls on either line
    if infinite == False and xp != None:
        if (A[0]-xp)*(xp-B[0]) < 0 or (C[0]-xp)*(xp-D[0]) < 0 or (A[1]-yp)*(yp-B[1]) < 0 or (C[1]-yp)*(yp-D[1]) < 0:
            xp,yp = None, None
    if xp != None:
        xp=float(xp)
    if yp != None:
        yp = float(yp)
    return (xp,yp)


def kvalue(p,a,b):
    """ 
    Determines the k-value of the point p on the line segment a-b. The k-value is the normalized
    location on the line ab, with k(a)=0 and k(b)=1

    Parameters
    ----------
    p : (x,y)
        Coordinates of a point. The point is assumed to be on the line through a & b.
    a : (x,y)
        Coordinates of the beginning of the line segment.
    b : (x,y)
        Coordinates of the end of the line segment.
    
    Returns
    ----------
    float
        The k-value of p on ab.
    """
    # The k-value can be computed from either x or y coordinates as long
    # as the values are not the same for a & b
    if(a==b):        # k-value is undefined
        return None
    elif a[0]==b[0]:    # use y-coordinates
        return (p[1]-a[1])/(b[1]-a[1])
    else: # use x-coordinates
        return (p[0]-a[0])/(b[0]-a[0])

    
def location(C,s,k):
    """
    Computes the location corresponding to the k-value along a segment of a polyline

    Parameters
    ----------
    C : [(x,y),...] list of tuples 
        The coordinates of the polyline.
    s : int
        The index of a segment on polyline C. Must be within [0,n-2]
    k : float
        The proportion from the start pt to the end pt of the segment.
    Returns
    ----------
    (x,y) : (float,float)
        The computed location.
    """
    x = C[s][0] + k*(C[s+1][0]-C[s][0])
    y = C[s][1] + k*(C[s+1][1]-C[s][1])
    return (x,y)    

    
def project_out(a1,a2,b1,b2):
    # returns k-value of perpendicular extension from b1 to a
    """ 
    Determines the k-value on [a1,a2] of the instersection between the line 
    through a1 and a2 and the line through b1 perpendicular to the line through
    b1 and b2. 

    Parameters
    ----------
    a1,a2,b1,b2 : (x,y)
        Coordinates of a point. 
    tryreverse: bool
        For internal use. Leave this alone.
    Returns
    ----------
    float
        The k-value of the projected out point on [a1,a2].
    """
        
    #Find the projection point on segment 'a' from a perpendicular line starting at b1
    H = project_pt_to_line(b1, a1, a2)
    # Find the k-value of the projection point
    HK = kvalue(H,a1,a2)
    #Find intersection of two segments
    I = intersection(a1, a2, b1, b2)
    # If there is no intersection, use k-value of projection point
    if I[0] == None:
        return HK
    else:
        #Find the distance from b1 to new point H to determine length
        H_B_Dist = distance(b1,H)
        #Find the distance from I to H to determine length
        H_I_Dist = distance(I,H)
        # if distance is zero, b1 is on a, so return its k-value
        if H_I_Dist == 0:
            return kvalue(H,a1,a2)
        else:
            #Find length of both unkown side of triangle
            K_H_Dist = (H_B_Dist**2) / H_I_Dist
            #Find the length of the entire 'a' segment
            L = distance(a1,a2)
            #Find the distance from intersection to the K point
            D = H_I_Dist + K_H_Dist
            #Find difference in k-values by dividing D by the overall length
            dk = D/L
            # Find the k-value of the intersection point
            IK = kvalue(I, a1, a2)
            # effective interval point is on far side of projection pt from intersection pt
            if HK < IK:
                k=IK-dk
            else:
                k = IK + dk
            return k

        
def rotate_pts(pts, origin, deg_cw):
    theta = __m.radians(deg_cw)
    cos_theta = __m.cos(-theta)
    sin_theta = __m.sin(-theta)
    # shift to origin
    p = origin
    s = [(x[0]-p[0],x[1]-p[1]) for x in pts]
    # rotate
    rotated_pts = []
    for p in s:
        x = p[0] * cos_theta - p[1] * sin_theta
        y = p[0] * sin_theta + p[1] * cos_theta
        rotated_pts.append((x,y))
    return rotated_pts


def rotate_horizontal(pts, pts2 = None):
    """
    Rotates two sets of points so the first and last point 
    of the first set form a horizontal line,
    with the first point on the left and at the origin.

    Parameters
    ----------
    pts : list of (x,y) tuples
        the points to be rotated.

    pts2 : list of (x,y) tuples (optional)
        a second set of points to be rotated

    Returns
    -------
    List of (x,y) tuples, or a tuple containing two lists of (x,y) tuples
        The rotated points.
    
    """
    # get point out from first point in horizontal direction
    hz_pt = (pts[0][0]+100,pts[0][1])
    # get angle between horizontal and line from first to last point
    a = angle(pts[0],pts[1],hz_pt)
    a = __m.degrees(a)
    print(a)
    # rotate
    r = rotate_pts(pts,pts[0], a+180)
    if pts2 != None:
        r2 = rotate_pts(pts2,pts[0], a+180)
    # shift to origin
    minx = min([x[0] for x in r])
    miny = min([x[1] for x in r])
    r = [(x[0]-minx,x[1]-miny) for x in r]
    if pts2 != None:
        r2 = [(x[0]-minx,x[1]-miny) for x in r2]
        return r,r2
    else:
        return r
