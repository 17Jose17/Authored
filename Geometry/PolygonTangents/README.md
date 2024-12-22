# Tangent Calculation Algorithm for Polygons

This algorithm calculates the tangents to a polygon given a set of external points. The time complexity of the algorithm is **O(n log n + k log n)**, where `n` is the number of vertices of the polygon, and `k` is the number of external points.

## Description

Given a polygon defined by `n` vertices and a set of `k` external points, the algorithm computes the tangents from each of the `k` points to the polygon. The tangents are printed as the output.

## Time Complexity

- **O(n log n + k log n)**

  - **O(n log n)**: This time complexity arises from the preprocessing step.
  - **O(k log n)**: To compute the tangents from the `k` external points.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
