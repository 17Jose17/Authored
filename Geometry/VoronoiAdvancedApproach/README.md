# **k-Order Voronoi Diagram Calculation**

This project implements an **incremental algorithm** to calculate the **k-order Voronoi diagram**. The algorithm works by dividing each region into the neighbors of the respective region.

## **Features**
- Efficient computation of k-order Voronoi diagrams using an incremental approach.
- Input: 
  - A set of points in 2D space.
  - An integer `k` representing the order of the Voronoi diagram.
- Output:
  - A 3D vector structure:
    1. **Index 0**: Indicates which input point owns the region (0-indexed).
    2. **Index 1**: Number of polygons in the region.
    3. **Index 2**: The polygons themselves, represented as collections of vertices.

---

## **License**

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for more details.

---

## **Contributing**

Contributions are welcome! Feel free to open issues or submit pull requests.

---

## **Acknowledgments**

This project is based on incremental algorithms for Voronoi diagram calculation and expands the method to k-order regions.
