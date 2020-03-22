import numpy as np
import math
from matplotlib import pyplot as plt


def split_medianX(pts_sortX, pts_sortY):
    ''' Splits the two arrays, one sorted by x-coords and one by y-coords, around the median x value
    
        If len(pts_sortX) is even: 
            the average of the two center elements is returned as the median x, 
            and the left and right arrays are of equal length.
        If len(pts_sortX) is odd: 
            the right subarray is one greater than the left
    '''
    n = len(pts_sortX)
    
    medianX = (pts_sortX[n//2 - 1][0] + pts_sortX[n//2][0]) / 2. if n%2==0 else pts_sortX[n//2][0]
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
    return math.sqrt(sum((pt2 - pt1)**2))
    

def pairwise_search(pts):
    ''' Brute-force O(n^2) closest-pair search among all points in the 
        passed array. Returns the closest pair and their Euclidean 
        distance.
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
    ''' Find the closest pair of points in pts_sortY_instrip that are
        less than max_dist of each other. If no such pair exists, returns
        None.
    '''
    Ys = pts_sortY_instrip[:,1]
    closest_dist = max_dist
    closest_pair = None

    # will run at most 6 times
    for i in range(len(pts_sortY_instrip)):
        j = i+1
        while (j < len(pts_sortY_instrip)) and (Ys[j] - Ys[i] < closest_dist):
            d_ij = distance(pts_sortY_instrip[i], pts_sortY_instrip[j])
            if d_ij < closest_dist:
                closest_dist = d_ij
                closest_pair = np.array([pts_sortY_instrip[i], pts_sortY_instrip[j]])
            j += 1

    return closest_pair


def plot_results(points, results, fig, ax, pause_len=1):
    for result in results:
        ax.clear()
        ax.set_xlim([0, 10])
        ax.set_ylim([0, 10])
        
        xmed = result['median']
        xminL = result['left min x']
        xmaxR = result['right max x']
        closest_pairL = result['left closest pair']
        closest_pairC = result['center closest pair']
        closest_pairR = result['right closest pair']

        ax.axvspan(xminL, xmed, alpha=0.05, color='red')
        ax.axvspan(xmed, xmaxR, alpha=0.05, color='blue')
        ax.axvline(xmed, ls='--', lw=1)
        ax.scatter(x=points[:,0], y=points[:,1], facecolors="None", edgecolors='k', lw=1, s=30)

        pairs = {}
        for i, (pair, c) in enumerate(zip([closest_pairL, closest_pairC, closest_pairR], ['r','g','b'])):
            if pair is not None:    
                x_vals = [pair[0][0], pair[1][0]]
                y_vals = [pair[0][1], pair[1][1]]
                pair_line = ax.plot(x_vals, y_vals, c=c, ls='-', lw=2)
                pair_Lpt = ax.scatter(x=x_vals[0], y=y_vals[0], facecolors=c, alpha=0.7, edgecolors="None", s=30)
                pair_Rpt = ax.scatter(x=x_vals[1], y=y_vals[1], facecolors=c, alpha=0.7, edgecolors="None", s=30)

                pair_dist = distance(pair[0], pair[1])
                pairs[pair_dist] = [pair_line, pair_Lpt, pair_Rpt]

        fig.canvas.draw()
        fig.canvas.flush_events()
        plt.pause(pause_len + 1)

        # remove all but closest pair from plot
        closest_pair = pairs.pop( min(pairs.keys()) )
        for key, val in pairs.items():
            ax.lines.remove(val[0][0])
            val[1].remove()
            val[2].remove()
        fig.canvas.draw()
        fig.canvas.flush_events()
        plt.pause(pause_len)

    # display overall closest pair
    ax.clear()
    ax.set_xlim([0,10])
    ax.set_ylim([0,10])
    ax.scatter(x=points[:,0], y=points[:,1], facecolors="None", edgecolors='k', lw=1, s=30)
    for pt in (closest_pair[1:]):
        pt.set_offset_position('data')
        x, y = pt.get_offsets()[0]
        ax.scatter(x=x, y=y, facecolors='y', edgecolors="None", alpha=0.8, s=30)
    #print(closest_pair[1].get_paths())
    #print(closest_pair[2].get_paths())
    fig.canvas.draw()
    fig.canvas.flush_events()

    plt.show(block=True)

#=============================================================================================

def _ClosestPair(pts_sortX, pts_sortY, results, level):
    ''' 
    Return values:
         (1) minimum distance among given points (infinity if <2 points)
         (2) points array sorted by y coordinate
    '''
    
    # base case
    if len(pts_sortX) <= 3:
        min_dist, closest_pair = pairwise_search(pts_sortX)
        return min_dist, closest_pair
    
    # recursively search left and right half-planes about median
    median, pts_sortX_L, pts_sortX_R, pts_sortY_L, pts_sortY_R \
        = split_medianX(pts_sortX, pts_sortY)
    
    delta_L, closest_pair_L = _ClosestPair(pts_sortX_L, pts_sortY_L, results, level+1)
    delta_R, closest_pair_R = _ClosestPair(pts_sortX_R, pts_sortY_R, results, level+1)
    
    # identify closest pair in left and right subsets
    if delta_L < delta_R:
        delta = delta_L
        closest_pair = closest_pair_L
    else:
        delta = delta_R
        closest_pair = closest_pair_R

    # merge left and right in order of y-coordinate
    #points_sortY = merge_sortY(pts_sortY_L, pts_sortY_R)
    
    # get all points within delta units of the median
    _, pts_sortY_instrip = pts_within_deltaX(pts_sortX, pts_sortY, delta, median)
    
    # check if any pair of points spanning the median is closer than delta units
    closest_pair_C = closest_pair_instrip(pts_sortY_instrip, delta)
    if closest_pair_C is not None:
        closest_pair = closest_pair_C
        delta = distance(closest_pair_C[0], closest_pair_C[1])

    results.append({'recursion level': level, 'median': median, 'left min x': pts_sortX_L[0][0],
        'left max x': pts_sortX_L[-1][0], 'left closest pair': closest_pair_L,
        'right min x': pts_sortX_R[0][0], 'right max x': pts_sortX_R[-1][0],
        'right closest pair': closest_pair_R, 'center closest pair': closest_pair_C})

    return delta, closest_pair
    
    
def ClosestPair(points):

    # presort by x and y coordinate
    pts_sortX = points[points[:,0].argsort()]
    pts_sortY = points[points[:,1].argsort()]
    
    results = []
    closest_dist, closest_pair = _ClosestPair(pts_sortX, pts_sortY, results, 1)

    return closest_dist, closest_pair, results

##############################################

if __name__ == "__main__":
    # def _clean_pairs(points):
    #     from collections import Counter

    #     points = np.asarray([p for p in points if None not in p])

    #     # prevent any duplicate x values with imperceptible jitter
    #     dup_Xs = {x: 0 for x, count in Counter(list(points[:, 0])).items() if count > 1}
    #     for i, x in enumerate(points[:, 0]):
    #         if x in dup_Xs:
    #             points[i, 0] += 1e-9 * dup_Xs[x]
    #             dup_Xs[x] += 1

    #     dup_Ys = {y: 0 for y, count in Counter(list(points[:, 1])).items() if count > 1}
    #     for i, y in enumerate(points[:, 1]):
    #         if y in dup_Ys:
    #             points[i, 1] += 1e-9 * dup_Ys[y]
    #             dup_Ys[y] += 1

    #     return points


    points1 = np.array(
        [[0.6322580645161291, 7.662337662337663],
         [1.6129032258064515, 6.987012987012987],
         [1.432258064516129, 5.012987012987013],
         [1.9225806451612906, 4.467532467532468],
         [3.470967741935484, 4.467532467532468],
         [7.341935483870968, 3.844155844155845],
         [8.141935483870967, 2.51948051948052],
         [9.380645161290323, 1.350649350649351],
         [8.890322580645162, 6.103896103896105],
         [6.154838709677419, 6.90909090909091],
         [8.529032258064516, 7.116883116883118],
         [6.903225806451612, 6.883116883116884],
         [2.980645161290323, 6.25974025974026],
         [5.199999999999999, 6.64935064935065],
         [5.225806451612904, 9.428571428571429],
         [3.341935483870968, 7.948051948051948],
         [4.3741935483870975, 7.0649350649350655],
         [4.348387096774193, 4.753246753246755],
         [5.251612903225807, 2.9870129870129873],
         [3.2903225806451615, 1.9740259740259742],
         [5.019354838709678, 1.9740259740259742],
         [5.638709677419355, 1.6363636363636367],
         [6.154838709677419, 1.5844155844155847],
         [5.974193548387097, 3.064935064935065],
         [2.232258064516129, 1.116883116883117],
         [1.1225806451612903, 2.12987012987013],
         [0.6064516129032258, 3.818181818181819],
         [1.5612903225806454, 3.636363636363637],
         [2.722580645161291, 3.064935064935065],
         [4.451612903225806, 0.9610389610389616],
         [6.516129032258064, 4.753246753246755]])

    points2 = np.array([[0,1], [2,1], [3,1], [5.5,1], [6,1], [7,1], [8,1], [10,1]])
    # [0,1], [2,1],  [3,1], [5.5,1]   |   [6,1], [7,1],   [8,1], [10,1]
    # [0,1], [2,1] | [3,1], [5.5,1]   |   [6,1], [7,1], | [8,1], [10,1]
    # 
    points = _clean_pairs(points1)
    closest_dist, closest_pair, results = ClosestPair(points)
    print(results)
    



