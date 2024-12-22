# Tangent Calculation Algorithm for Polygons

This algorithm calculates the tangents to a polygon given a set of external points. The time complexity of the algorithm is **O(n log n + k log n)**, where `n` is the number of vertices of the polygon, and `k` is the number of external points.

## Description

Given a polygon defined by `n` vertices and a set of `k` external points, the algorithm computes the tangents from each of the `k` points to the polygon. The tangents are printed as the output.

## Time Complexity

- **O(n log n + k log n)**

  - **O(n log n)**: To sort the vertices of the polygon.
  - **O(k log n)**: To compute the tangents from the `k` external points.

## How to Use

1. The program first reads the number of vertices `n` of the polygon.
2. Then, it reads the `n` vertices of the polygon, represented by their (x, y) coordinates.
3. Next, it reads the number `k` of external points.
4. Then, it reads the coordinates of the `k` external points.
5. Finally, the program prints the tangents from each of the `k` points to the polygon.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
