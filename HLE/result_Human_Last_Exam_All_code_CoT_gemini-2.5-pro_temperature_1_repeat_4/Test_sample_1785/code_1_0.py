import math

def calculate_rectangle_edges(corner1, corner2):
    """
    Calculates the number of edges (perimeter) of a lattice rectangle
    aligned with the coordinate axes.
    The rectangle is defined by two opposite corners.
    """
    # The number of edges is the perimeter of the rectangle.
    # For a lattice rectangle, this is 2 * (width + height).
    
    # Calculate the sum of side lengths along each dimension.
    # For a rectangle, only two dimensions will have a non-zero length.
    length_sum = abs(corner1[0] - corner2[0]) + \
                 abs(corner1[1] - corner2[1]) + \
                 abs(corner1[2] - corner2[2])
                 
    # The perimeter is twice the sum of the two side lengths.
    return 2 * length_sum

# The minimal known configuration for Borromean rings on a cubic lattice
# consists of three mutually perpendicular rectangles.

# Component 1 (K1): A rectangle in the xy-plane (z=0)
c1_corner1 = (1, 2, 0)
c1_corner2 = (-1, -2, 0)
edges_k1 = calculate_rectangle_edges(c1_corner1, c1_corner2)

# Component 2 (K2): A rectangle in the yz-plane (x=0)
c2_corner1 = (0, 1, 2)
c2_corner2 = (0, -1, -2)
edges_k2 = calculate_rectangle_edges(c2_corner1, c2_corner2)

# Component 3 (K3): A rectangle in the xz-plane (y=0)
c3_corner1 = (2, 0, 1)
c3_corner2 = (-2, 0, -1)
edges_k3 = calculate_rectangle_edges(c3_corner1, c3_corner2)

# Calculate the total number of edges
total_edges = edges_k1 + edges_k2 + edges_k3

# Print the result as an equation
print("The minimum total number of edges is achieved by three components with the following number of edges:")
print(f"{edges_k1} + {edges_k2} + {edges_k3} = {total_edges}")
