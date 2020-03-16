import numpy as np
import math
from matplotlib import pyplot as plt
import time


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


def plot_results(points, results, fig, ax, pause_len=1.5):
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
        ax.scatter(x=x, y=y, facecolors='y')
    #print(closest_pair[1].get_paths())
    #print(closest_pair[2].get_paths())
    fig.canvas.draw()
    fig.canvas.flush_events()

    plt.show(block=True)

#=============================================================================================

def _ClosestPair(pts_sortX, pts_sortY, results):
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
    
    delta_L, closest_pair_L = _ClosestPair(pts_sortX_L, pts_sortY_L, results)
    delta_R, closest_pair_R = _ClosestPair(pts_sortX_R, pts_sortY_R, results)
    
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

    results.append({'median': median, 'left min x': pts_sortX_L[0][0],
        'left max x': pts_sortX_L[-1][0], 'left closest pair': closest_pair_L,
        'right min x': pts_sortX_R[0][0], 'right max x': pts_sortX_R[-1][0],
        'right closest pair': closest_pair_R, 'center closest pair': closest_pair_C})

    return delta, closest_pair
    
    
def ClosestPair(points):

    # presort by x and y coordinate
    pts_sortX = points[points[:,0].argsort()]
    pts_sortY = points[points[:,1].argsort()]
    
    results = []
    closest_dist, closest_pair = _ClosestPair(pts_sortX, pts_sortY, results)

    return closest_dist, closest_pair, results

##############################################

if __name__ == "__main__":
    points1 = np.array(
        [[0.504032258064516, 0.38419913419913465],
         [0.786290322580645, 3.360389610389611],
         [1.0181451612903225, 3.37391774891775],
         [1.411290322580645, 3.8879870129870135],
         [7.006048387096774, 7.621753246753247],
         [5.120967741935484, 6.5530303030303045],
         [2.7217741935483875, 7.824675324675326],
         [8.104838709677418, 5.267857142857144],
         [8.316532258064516, 2.7245670995671003],
         [3.4274193548387095, 3.37391774891775],
         [6.350806451612904, 4.686147186147187],
         [7.772177419354838, 6.201298701298702],
         [8.689516129032258, 7.337662337662339],
         [8.16532258064516, 6.174242424242426],
         [4.818548387096774, 7.851731601731602],
         [5.181451612903226, 6.3365800865800885],
         [4.485887096774194, 4.983766233766234],
         [5.473790322580646, 3.7256493506493515],
         [5.01008064516129, 1.4799783549783556],
         [0.282258064516129, 9.691558441558442]])

    points2 = np.array([[0,1], [2,1], [4.5,1], [5.5,1], [8,1], [10,1]])

    closest_dist, closest_pair, results = ClosestPair(points1)

    fig = plt.figure(figsize=(10,10))
    plt.ion()
    ax = fig.add_subplot(111)
    ax.set_xlim([-1, 11])
    ax.set_ylim([-1, 11])
    plot_results(points1, results, fig, ax)

