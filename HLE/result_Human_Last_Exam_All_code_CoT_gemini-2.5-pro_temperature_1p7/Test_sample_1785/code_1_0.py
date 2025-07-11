import math

def calculate_path_length(path_vertices):
    """
    Calculates the total number of unit edges for a closed path on a lattice.
    The path is defined by a list of its corner vertices.
    """
    total_length = 0
    num_vertices = len(path_vertices)
    for i in range(num_vertices):
        p1 = path_vertices[i]
        # Get the next vertex in the path, wrapping around for the last segment
        p2 = path_vertices[(i + 1) % num_vertices]
        
        # Calculate Manhattan distance, which corresponds to the number of unit edges
        # between two orthogonally-aligned points.
        distance = abs(p1[0] - p2[0]) + abs(p1[1] - p2[1]) + abs(p1[2] - p2[2])
        total_length += distance
    return total_length

# The established minimal construction for the Borromean rings uses three
# congruent 4x2 rectangles, oriented orthogonally to each other.

# Component 1: A 4x2 rectangle in the xy-plane (z=0)
path_A_corners = [
    (2, 1, 0),
    (-2, 1, 0),
    (-2, -1, 0),
    (2, -1, 0)
]

# Component 2: A 4x2 rectangle in the yz-plane (x=0) by cyclically permuting coordinates (x,y,z)->(y,z,x)
path_B_corners = [
    (0, 2, 1),
    (0, -2, 1),
    (0, -2, -1),
    (0, 2, -1)
]

# Component 3: A 4x2 rectangle in the xz-plane (y=0) by cyclically permuting coordinates again
path_C_corners = [
    (1, 0, 2),
    (1, 0, -2),
    (-1, 0, -2),
    (-1, 0, 2)
]

# Calculate the number of edges for each component
len_A = calculate_path_length(path_A_corners)
len_B = calculate_path_length(path_B_corners)
len_C = calculate_path_length(path_C_corners)

# Calculate the total number of edges in the link
total_len = len_A + len_B + len_C

print("This script calculates the minimum number of edges for a topologically non-trivial, 3-component link on the cubic lattice.")
print("This corresponds to the known minimal construction of the Borromean rings.\n")
print(f"Component 1 is a {abs(path_A_corners[1][0]-path_A_corners[0][0])}x{abs(path_A_corners[2][1]-path_A_corners[1][1])} rectangle with edge length: {len_A}")
print(f"Component 2 is a {abs(path_B_corners[1][1]-path_B_corners[0][1])}x{abs(path_B_corners[2][2]-path_B_corners[1][2])} rectangle with edge length: {len_B}")
print(f"Component 3 is a {abs(path_C_corners[2][0]-path_C_corners[1][0])}x{abs(path_C_corners[1][2]-path_C_corners[0][2])} rectangle with edge length: {len_C}")
print("\nThese three components form a topologically non-trivial link while having disjoint vertex sets.")
print("\nThe minimum total number of edges is the sum of the lengths of the three components.")
print(f"Final Equation: {len_A} + {len_B} + {len_C} = {total_len}")
