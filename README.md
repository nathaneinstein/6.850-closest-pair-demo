## 6.850 Computational Geometry: Problem Set 1: *closest pair implementation*
Demo of the divide-and-conquer algorithm for finding the closest pair of points in the Euclidean plane.

## Usage and Dependencies

- NumPy
- Matplotlib

To run:
```bash
./python main.py
```

## Pseudocode

Store points in sorted order by both x- and y-coordinate. 

*Recursively...*
 1. Divide the points about the median x-coordinate into 2 approximately equal-sized subsets. (This operation involves checking each point in the array that stores the points in order of their y-coordinate values to avoid re-sorting them at each level of the recursion tree).
 2. Search for the closest pair in the left and right subsets, and find the minimum ($\delta$) of the distance between these two candidate closest pairs
 3. Check for a potentially closer pair of two points divided between the two subsets (i.e. one to the left of the median, one to the right). Only those points that fall within $\delta$ units of the median line need to be checked, and for each of these points, only 5 vertically adjacent neighbors need to be considered.

The base case occurs when there are 3 or fewer points in the subset, in which case the closest pair can be found through a brute-force comparison of all pairs.

---

The demo visualizes the "candidate" closest pairs at each level in the recursion tree. 

The left and right subsets at each level are highlighted with red and blue backgrounds, respectively, and the closest pairs identified in the left and right subtrees are connected with red and blue lines, respectively. If a closer pair of points spans the median line, that pair is also displayed with a green connecting line. After stepping through each level of the recursion tree, the closest pair overall is highlighted in yellow. 

An example is shown below:



![Demo](https://github.com/nathaneinstein/6.850-closest-pair-demo/blob/master/algdemo.gif)
