# Polyhedron Cutting Algorithm

### Description
This algorithm is designed to cut a convex or non-convex polyhedron with any arbitrary plane. It operates with linear complexity, making it highly efficient for computational geometry tasks. The algorithm takes as input:

- A polyhedron defined by its number of faces and, for each face, the number of vertices and their coordinates.
- A plane defined by its normal vector and the value `d` in the plane equation.

### Key Features:
- Handles both convex and non-convex polyhedra.
- Outputs the resulting sub-polyhedra or modified geometry after the cut.
- Optimized for linear time complexity.

This implementation is suitable for applications such as 3D modeling, physics simulations, and geometric computations.

### Usage
To use the algorithm:
1. Define the polyhedron by specifying:
   - The total number of faces.
   - For each face, the number of vertices and their coordinates in 3D space.
2. Specify the cutting plane using:
   - A normal vector `(a, b, c)`.
   - The value `d` in the plane equation `ax + by + cz + d = 0`.
3. Run the algorithm to obtain the resulting geometry after the cut.

Refer to the source code for detailed examples and documentation on input and output formats.

### License
This project is licensed under the **Creative Commons Attribution-NonCommercial-NoDerivatives 3.0 Unported License (CC BY-NC-ND 3.0)**.

#### License Summary:
1. **Attribution**: You must give appropriate credit to the original author.
2. **Non-Commercial**: You may not use this code for commercial purposes.
3. **No Derivatives**: You may not distribute modified versions of this code.

For the full license text, see the `LICENSE` file or visit: [https://creativecommons.org/licenses/by-nc-nd/3.0/](https://creativecommons.org/licenses/by-nc-nd/3.0/)

**Author**: [José Ehécatl Mejía Yáñez]  
**Date**: December 21, 2024
