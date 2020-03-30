import tkinter
from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg, NavigationToolbar2Tk)
from matplotlib.backend_bases import key_press_handler
from matplotlib.figure import Figure
import numpy as np
from collections import Counter

import closest_pair as cp


root = tkinter.Tk()
root.wm_title("Closest-Pair Algorithm Demo")

fig = Figure(figsize=(5, 5), dpi=100)
ax = fig.add_subplot(111)
ax.set_xlim([0, 10])
ax.set_ylim([0, 10])

canvas = FigureCanvasTkAgg(fig, master=root)
canvas.draw()
canvas.get_tk_widget().pack(side=tkinter.TOP, fill=tkinter.BOTH, expand=1)

toolbar = NavigationToolbar2Tk(canvas, root)
toolbar.update()
canvas.get_tk_widget().pack(side=tkinter.TOP, fill=tkinter.BOTH, expand=1)


#---------------------------------------------------
# Allow user to specify points by clicking on canvas
#---------------------------------------------------
points = []
def _add_pt(event):
    ax.scatter(event.xdata, event.ydata, facecolors='none', edgecolors='k', s=30)
    fig.canvas.draw()
    points.append([event.xdata, event.ydata])

canvas.mpl_connect('button_press_event', _add_pt)


#-----------------------------------------------------
# Add button to run closest-pair alg on points entered
#-----------------------------------------------------
def _clean_pairs(points):
	# remove [None, None] pairs from clicking outside of plot bounds
	points = np.asarray([p for p in points if None not in p])
	
	# ensure no duplicate x values by adding imperceptible jitter
	dup_Xs = {x:0 for x, count in Counter(list(points[:,0])).items() if count > 1}
	for i, x in enumerate(points[:,0]):
		if x in dup_Xs:
			points[i,0] += 1e-9 * dup_Xs[x]
			dup_Xs[x] += 1

	# ensure no duplicate y values by adding imperceptible jitter
	dup_Ys = {y:0 for y, count in Counter(list(points[:,1])).items() if count > 1}
	for i, y in enumerate(points[:,1]):
		if y in dup_Ys:
			points[i,1] += 1e-9 * dup_Ys[y]
			dup_Ys[y] += 1
	
	return points

def _find_closest_pair(points, fig, ax):   
    print("\nPoints entered:\n---------------\n", points)
    points = _clean_pairs(points)
    closest_dist, closest_pair, results = cp.ClosestPair(points)
    with np.printoptions(precision=3):
	    print("\nClosest-pair results:\n---------------------\n", results)
    cp.plot_results(points, results, fig, ax)

run_alg_button = tkinter.Button(master=root, text="Find closest pair", 
                                command=lambda: _find_closest_pair(points, fig, ax)) 
run_alg_button.pack(side=tkinter.BOTTOM)


tkinter.mainloop()