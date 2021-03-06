import numpy as np
import math
from matplotlib import pyplot as plt


def split_medianX(pts_sortX, pts_sortY):
    ''' Splits the two arrays, one sorted by x-coords and one by y-coords, 
        around the median x value
    
        If len(pts_sortX) is even: 
            the average of the two center elements is returned as the median x, 
            and the left and right arrays are of equal length.
        If len(pts_sortX) is odd: 
            the right subarray is one greater than the left
    '''
    n = len(pts_sortX)
    
    if n%2==0:
        medianX = (pts_sortX[n//2 - 1][0] + pts_sortX[n//2][0]) / 2. 
    else:
        medianX = pts_sortX[n//2][0]
    pts_sortX_L = pts_sortX[:n//2] 
    pts_sortX_R = pts_sortX[n//2:]
    
    pts_sortY_L = np.empty(pts_sortX_L.shape)
    pts_sortY_R = np.empty(pts_sortX_R.shape)
    ptrL = 0; ptrR = 0
    for pt in pts_sortY:
        if pt[0] < medianX:
            pts_sortY_L[ptrL] = pt
            ptrL += 1
        else:
            pts_sortY_R[ptrR] = pt
            ptrR += 1
    
    return medianX, pts_sortX_L, pts_sortX_R, pts_sortY_L, pts_sortY_R


def distance(pt1, pt2): 
    ''' Return the Euclidean distance between points pt1 and pt2 '''
    return math.sqrt(sum((pt2 - pt1)**2))
    

def pairwise_search(pts):
    ''' Performs a brute-force O(n^2) search for closest pair among
        all points in the passed array. Returns the closest pair
        and their Euclidean distance.
    '''
    min_dist = math.inf
    closest_pair = None

    n = len(pts)
    for i in range(n):
        for j in range(i+1,n):
            d = distance(pts[i], pts[j])
            if d < min_dist:
                min_dist = d
                closest_pair = np.array([pts[i], pts[j]])

    return min_dist, closest_pair


def pts_within_deltaX(pts_sortX, pts_sortY, delta, x):
    ''' Returns points from pts_sortX and pts_sortY that fall within 
        delta distance of the specified x value.
    '''
    # take advantage of sorted Xs to find strip edges in O(logn)
    lower = np.searchsorted(pts_sortX[:,0], x - delta - 1e-6) # add eps to include pt at x-delta
    upper = np.searchsorted(pts_sortX[:,0], x + delta, side="right")
    pts_sortX_instrip = pts_sortX[lower:upper]
    
    pts_sortY_instrip = np.empty(pts_sortX_instrip.shape)
    i = 0
    for pt in pts_sortY:
        if (pt[0] >= x-delta) and (pt[0] <= x+delta):
            pts_sortY_instrip[i] = pt
            i += 1
            
    return pts_sortX_instrip, pts_sortY_instrip
    
    
def closest_pair_instrip(pts_sortY_instrip, max_dist):
    ''' Find the closest pair of points in pts_sortY_instrip -- the
        subset of all points that fall within max_dist of the median line,
        sorted by y coordinate --  that are < max_dist of each other. 
        If no such pair exists, returns None.
    '''
    Ys = pts_sortY_instrip[:,1]
    closest_dist = max_dist
    closest_pair = None
    pairs_checked = []

    # inner loop will run at most 5 times. This fact is key to demonstrating 
    # the algorithm runs in O(nlogn) rather than O(n^2)
    for i in range(len(pts_sortY_instrip)):
        j = i+1
        while (j < len(pts_sortY_instrip)) and (Ys[j] - Ys[i] < closest_dist):
            pt_i = pts_sortY_instrip[i]
            pt_j = pts_sortY_instrip[j]
            pairs_checked.append((pt_i, pt_j))
            d_ij = distance(pt_i, pt_j)
            if d_ij < closest_dist:
                closest_dist = d_ij
                closest_pair = np.array([pt_i, pt_j])
            j += 1

    return closest_pair, pairs_checked


def plot_results(points, results, fig, ax, pause_len=2):

    def refresh_plot(pause=pause_len):
        fig.canvas.draw()
        fig.canvas.flush_events()
        plt.pause(pause)


    for result in results:
        ax.clear()
        ax.set_xlim([0, 10])
        ax.set_ylim([0, 10])
        
        xmed = result['median']
        xminL = result['L min x']
        xmaxR = result['R max x']
        closest_pairL = result['L closest pair']
        closest_pairC = result['center closest pair']
        closest_pairR = result['R closest pair']
        pairs_checked = result['center pairs checked']

        ax.axvspan(xminL, xmed, alpha=0.05, color='red')
        ax.axvspan(xmed, xmaxR, alpha=0.05, color='blue')
        ax.axvline(xmed, ls='--', lw=2)
        ax.scatter(x=points[:,0], y=points[:,1], facecolors="None", edgecolors='k', lw=1, s=30)

        # draw closest pairs in left and right subsets
        pairs = {}
        for i, (pair, c) in enumerate(zip([closest_pairL, closest_pairR], ['r','b'])):
            x_vals = [pair[0][0], pair[1][0]]
            y_vals = [pair[0][1], pair[1][1]]
            pair_line = ax.plot(x_vals, y_vals, c=c, ls=':', lw=1)
            pair_Lpt = ax.scatter(x=x_vals[0], y=y_vals[0], facecolors=c, alpha=0.7, edgecolors="None", s=30)
            pair_Rpt = ax.scatter(x=x_vals[1], y=y_vals[1], facecolors=c, alpha=0.7, edgecolors="None", s=30)

            pair_dist = distance(pair[0], pair[1])
            pairs[pair_dist] = [pair_line, pair_Lpt, pair_Rpt]
        refresh_plot()

        # remove all but closest pair over left and right subsets
        delta = min(pairs.keys())
        closest_pair = pairs.pop(delta)
        for key, val in pairs.items():
            ax.lines.remove(val[0][0])
            val[1].remove()
            val[2].remove()
        refresh_plot()

        # show pair comparisons within 2*delta strip about median
        ax.axvspan(xmed - delta, xmed + delta, alpha=0.1, color='green')
        refresh_plot()
        closest_pair_ctr = None
        for pair in pairs_checked:
            x_vals = [pair[0][0], pair[1][0]]
            y_vals = [pair[0][1], pair[1][1]]
            pair_Lpt = ax.scatter(x=x_vals[0], y=y_vals[0], facecolors="None", linewidths=2, edgecolors='k', s=30)
            pair_Rpt = ax.scatter(x=x_vals[1], y=y_vals[1], facecolors="None", linewidths=2, edgecolors='k', s=30)
            refresh_plot()
            pair_Lpt.remove()
            pair_Rpt.remove()
            refresh_plot()

        # if closest pair is in strip about median, show it
        if closest_pairC is not None:
            ax.lines.remove(closest_pair[0][0]); closest_pair[1].remove(); closest_pair[2].remove()
            x_vals = [closest_pairC[0][0], closest_pairC[1][0]]
            y_vals = [closest_pairC[0][1], closest_pairC[1][1]]
            pair_line = ax.plot(x_vals, y_vals, c='g', ls=':', lw=1)
            pair_Lpt = ax.scatter(x=x_vals[0], y=y_vals[0], facecolors='g', alpha=0.7, edgecolors="None", s=30)
            pair_Rpt = ax.scatter(x=x_vals[1], y=y_vals[1], facecolors='g', alpha=0.7, edgecolors="None", s=30)
            closest_pair = [pair_line, pair_Lpt, pair_Rpt]
            refresh_plot()

    # display overall closest pair
    ax.clear()
    ax.set_xlim([0,10])
    ax.set_ylim([0,10])
    ax.scatter(x=points[:,0], y=points[:,1], facecolors="None", edgecolors='k', lw=1, s=30)
    for pt in (closest_pair[1:]):
        pt.set_offset_position('data')
        x, y = pt.get_offsets()[0]
        ax.scatter(x=x, y=y, facecolors='y', edgecolors="k", s=45, linewidths=2)
    refresh_plot(0)

    plt.show(block=True)

#===============================================================================

def _ClosestPair(pts_sortX, pts_sortY, results, level):
    ''' 
    Helper function for recursive ClosestPair function
    
    Parameters:
        pts_sortX:  points sorted by x-coordinate
        pts_sortY:  points (same set as in pts_sortX) sorted by y-coordinate
        results:    list of dictionaries storing the results at each level in 
                    the recursion tree.
        level:      current level in the recursion tree

    Return values:
        1. Minimum distance of all pairs of points
        2. Coordinates of the closest pair of points
    '''
    
    # base case
    if len(pts_sortX) <= 3:
        return pairwise_search(pts_sortX)
    
    # split points about median x value and recursively search these subsets
    median, pts_sortX_L, pts_sortX_R, pts_sortY_L, pts_sortY_R \
        = split_medianX(pts_sortX, pts_sortY)
    
    delta_L, closest_pair_L = _ClosestPair(pts_sortX_L, pts_sortY_L, results, level+1)
    delta_R, closest_pair_R = _ClosestPair(pts_sortX_R, pts_sortY_R, results, level+1)
    
    # identify which of left- or right-subset candidate closest pairs is closer
    if delta_L < delta_R:
        delta = delta_L
        closest_pair = closest_pair_L
    else:
        delta = delta_R
        closest_pair = closest_pair_R

    # get all points within delta-unit strip about the median
    _, pts_sortY_instrip = pts_within_deltaX(pts_sortX, pts_sortY, delta, median)
    
    # check if any pair of points spanning the median is closer than delta units
    closest_pair_C, pairs_checked = closest_pair_instrip(pts_sortY_instrip, delta)
    if closest_pair_C is not None:
        closest_pair = closest_pair_C
        delta = distance(closest_pair_C[0], closest_pair_C[1])

    results.append({'recursion level': level, 'median': median,
        'L min x': pts_sortX_L[0][0], 'L max x': pts_sortX_L[-1][0], 
        'L closest pair': closest_pair_L, 'R min x': pts_sortX_R[0][0], 
        'R max x': pts_sortX_R[-1][0], 'R closest pair': closest_pair_R, 
        'center closest pair': closest_pair_C, 
        'center pairs checked': pairs_checked})

    return delta, closest_pair
    
    
def ClosestPair(points):
    ''' Run the divide-and-conquer closest-pair algorithm on the set of passed
        points in a 2D plane. 
    '''

    # presort by x and y coordinate
    pts_sortX = points[points[:,0].argsort()]
    pts_sortY = points[points[:,1].argsort()]
    
    results = []
    closest_dist, closest_pair = _ClosestPair(pts_sortX, pts_sortY, results, 1)

    return closest_dist, closest_pair, results



