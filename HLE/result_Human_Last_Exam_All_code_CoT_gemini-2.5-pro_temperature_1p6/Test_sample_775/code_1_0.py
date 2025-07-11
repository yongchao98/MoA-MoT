import numpy as np

def demonstrate_components(grid_size=15, n_components=4):
    """
    Demonstrates a construction on a discrete grid where the intersection
    of two connected sets has n_components.

    - Grid legend:
      'A': Belongs to set A only
      'B': Belongs to set B only
      'I': Belongs to the Intersection (A and B)
      '.': Empty
    """
    if n_components < 1:
        print("Number of components must be at least 1.")
        return

    # Initialize grid
    grid = np.full((grid_size, grid_size), 'B', dtype=str)

    # 1. Define the intersection set K (n_components points)
    k_points = []
    for i in range(n_components):
        # Place points in a zig-zag pattern
        r = (grid_size // (n_components + 1)) * (i + 1)
        c = (grid_size // 4) * (1 + 2 * (i % 2))
        if 0 < r < grid_size and 0 < c < grid_size:
            k_points.append((r, c))
    
    # Recalculate n_components in case some points were out of bounds
    n_components = len(k_points)
    if n_components == 0:
        print(f"Grid size {grid_size} is too small to place components.")
        return

    # 2. Define A as a path connecting the points in K
    # A will be the path plus the intersection points.
    # A_prime is the path without the intersection points.
    
    # Sort points to make a clear path
    k_points.sort()
    
    for i in range(len(k_points) - 1):
        p1 = k_points[i]
        p2 = k_points[i+1]
        
        # Draw a simple L-shaped path from p1 to p2
        r1, c1 = p1
        r2, c2 = p2
        
        # Horizontal part
        for c in range(min(c1, c2), max(c1, c2) + 1):
            grid[r1, c] = 'A'
        # Vertical part
        for r in range(min(r1, r2), max(r1, r2) + 1):
            grid[r, c2] = 'A'

    # Mark the intersection points
    for r, c in k_points:
        grid[r, c] = 'I'
        
    # Set A is all 'A' and 'I' cells. It is connected by construction.
    # Set B is all 'B' and 'I' cells. It is the background, which is also connected.
    # The intersection is the set of 'I' cells.
    
    print(f"Construction for a grid where the intersection has {n_components} components:\n")
    for row in grid:
        print(" ".join(row))
    
    print(f"\nFinal Equation: In this discrete analogy, the number of components is {n_components}.")
    print("In the continuous plane, this number can be arbitrarily large.")


demonstrate_components(grid_size=17, n_components=4)
