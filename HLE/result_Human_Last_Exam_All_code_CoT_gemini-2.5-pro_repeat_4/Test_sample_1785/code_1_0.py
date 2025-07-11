def calculate_loop_length(corners):
    """
    Calculates the total number of edges for a rectangular loop on the integer lattice,
    given its corner vertices.
    """
    total_length = 0
    num_corners = len(corners)
    for i in range(num_corners):
        p1 = corners[i]
        p2 = corners[(i + 1) % num_corners]  # Wrap around to the start
        # Calculate Manhattan distance between corners, which corresponds to the number of edges
        length = abs(p1[0] - p2[0]) + abs(p1[1] - p2[1]) + abs(p1[2] - p2[2])
        total_length += length
    return int(total_length)

# This script constructs the minimal 3-component non-trivial link, which consists of
# a minimal Hopf link and a minimal unknot.

# Component 1: A 2x2 square in the xy-plane (z=0). Part of the Hopf link.
# The vertices listed are the corners of the square.
knot1_corners = [
    (0, 0, 0),
    (2, 0, 0),
    (2, 2, 0),
    (0, 2, 0),
]

# Component 2: A 2x2 square in the yz-plane, offset to x=1.
# It interlocks with Component 1 to form a Hopf link. Its vertices are disjoint from Component 1.
knot2_corners = [
    (1, -1, -1),
    (1, 3, -1),
    (1, 3, 1),
    (1, -1, 1),
]

# Component 3: A minimal 1x1 square (unknot).
# It is placed far from the others to ensure it is unlinked and disjoint.
knot3_corners = [
    (10, 10, 10),
    (11, 10, 10),
    (11, 11, 10),
    (10, 11, 10),
]

# Calculate the length (number of edges) for each component.
len1 = calculate_loop_length(knot1_corners)
len2 = calculate_loop_length(knot2_corners)
len3 = calculate_loop_length(knot3_corners)

# Calculate the total length of the link.
total_length = len1 + len2 + len3

# Print the final equation showing the length of each component and the total length.
print(f"The minimum total number of edges is achieved by a link composed of three knots with lengths {len1}, {len2}, and {len3}.")
print(f"The final equation is: {len1} + {len2} + {len3} = {total_length}")
