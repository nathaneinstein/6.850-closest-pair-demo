import tkinter
from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg, NavigationToolbar2Tk)
# Implement the default Matplotlib key bindings.
from matplotlib.backend_bases import key_press_handler
from matplotlib.figure import Figure
import numpy as np

import closest_pair as cp


root = tkinter.Tk()
root.wm_title("Closest-Pair Algorithm Demo")

fig = Figure(figsize=(5, 5), dpi=100)
ax = fig.add_subplot(111)
ax.set_xlim([0, 10])
ax.set_ylim([0, 10])


canvas = FigureCanvasTkAgg(fig, master=root)  # A tk.DrawingArea.
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
#    print('button=%d, x=%d, y=%d, xdata=%f, ydata=%f' %
#          (event.button, event.x, event.y, event.xdata, event.ydata))
    ax.scatter(event.xdata, event.ydata, facecolors='none', edgecolors='k', s=30)
    fig.canvas.draw()
    points.append([event.xdata, event.ydata])

canvas.mpl_connect('button_press_event', _add_pt)

#---------------------------------------------------
# Button to trigger closest-pair algorithm
#---------------------------------------------------

def _clean_pairs(points):
	from collections import Counter
	# remove [None, None] pairs from clicking outside of plot bounds
	points = np.asarray([p for p in points if None not in p])
	
	# prevent any duplicate x values with imperceptible jitter
	dup_Xs = {x:0 for x, count in Counter(list(points[:,0])).items() if count > 1}
	for i, x in enumerate(points[:,0]):
		if x in dup_Xs:
			points[i,0] += 1e-16 * dup_Xs[x]
			dup_Xs[x] += 1
	
	dup_Ys = {y:0 for y, count in Counter(list(points[:,1])).items() if count > 1}
	for i, y in enumerate(points[:,1]):
		if y in dup_Ys:
			points[i,1] += 1e-16 * dup_Ys[y]
			dup_Ys[y] += 1
	
	return points

def _find_closest_pair(points, fig, ax):
    print(points)
    points = _clean_pairs(points)
    closest_dist, closest_pair, results = cp.ClosestPair(points)
    print(results)
    cp.plot_results(points, results, fig, ax)


run_alg_button = tkinter.Button(master=root, text="Find closest pair", 
                                command=lambda: _find_closest_pair(points, fig, ax)) 
run_alg_button.pack(side=tkinter.BOTTOM)


#---------------------------------------------------
# Quit button
#---------------------------------------------------
# def _quit():
#     root.quit()
#     root.destroy()                  

# quit_button = tkinter.Button(master=root, text="Quit", command=_quit)
# quit_button.pack(side=tkinter.BOTTOM)

tkinter.mainloop()