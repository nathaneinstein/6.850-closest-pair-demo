## **6.850 Computational Geometry**, Problem Set 1: *closest pair implementation*
Demo of the divide-and-conquer algorithm for finding the closest pair of points in the Euclidean plane.

## Usage and Dependencies

Python dependencies:
- NumPy
- Matplotlib
- Python 3.* (tested with version 3.7)

To run:
```bash
python main.py
```

## Pseudocode

Store points in sorted order by both x- and y-coordinate. 

*Recursively...*
 1. Divide the points about the median x-coordinate into 2 approximately equal-sized subsets. (This operation involves checking each point in the array that stores the points in order of their y-coordinate values to avoid re-sorting them at each level of the recursion tree).
 2. Search for the closest pair in the left and right subsets, and find the minimum (𝛿) of the distance between these two candidate closest pairs
 3. Check for a potentially closer pair of two points divided between the two subsets (i.e. one to the left of the median, one to the right). Only those points that fall within 𝛿 units of the median line need to be checked, and for each of these points, only 5 vertically adjacent neighbors need to be considered.

The base case occurs when there are 3 or fewer points in the subset, in which case the closest pair can be found through a brute-force comparison of all pairs.

---

The demo visualizes the "candidate" closest pairs at each level in the recursion tree. In particular, for each subtree, the range of x-values in the left and right subsets about the median x-coordinate are highlighted with red and blue backgrounds, respectively. The median x-value is illustrated with a vertical dotted line dividing these two subregions. The demo:

1. Displays the closest pair found recursively for each of these two subregions by highlighting the points in red or blue (corresponding to which subtree they fall in) and connecting each pair with a dotted line. 
2. The further of these two pairs is then deselected, and a green strip illustrating the range of x-values ±𝛿 units of the median is displayed.
3. The demo then visualizes the pairwise comparisons that are carried out for all points in this strip whose y-coordinates differ by at most 𝛿 by sequentially bolding each pair of points. If the Euclidean distance between one of these pairs is found to be less than 𝛿, the new closest pair is visualized in green. 

After stepping through each level of the recursion tree, the closest pair overall is highlighted in yellow. 

An example is shown below:



![Demo](https://github.com/nathaneinstein/6.850-closest-pair-demo/blob/master/algdemo.gif)
