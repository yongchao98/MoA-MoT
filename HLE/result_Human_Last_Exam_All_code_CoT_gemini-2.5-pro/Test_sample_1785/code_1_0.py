def calculate_path_length(path):
    """Calculates the number of edges in a path where the last point connects to the first."""
    # The number of edges in a closed loop is the number of vertices listed
    # assuming the path is given as p1 -> p2 -> ... -> pN -> p1.
    # Our path list is [p1, p2, ..., pN, p1], so it has N+1 points for N edges.
    return len(path) - 1

def main():
    """
    This program demonstrates the minimal edge construction for a 3-component
    nontrivial link (Borromean rings) on a 3D integer lattice.
    The construction uses three 6-edge loops.
    """

    # The coordinates for the three disjoint 6-edge loops (knots).
    # Each path is defined as a list of vertices, closing the loop at the end.
    
    # Knot 1
    knot1_path = [
        (0, 0, 1), (1, 0, 1), (1, 0, 0),
        (1, 1, 0), (0, 1, 0), (0, 0, 0),
        (0, 0, 1)
    ]

    # Knot 2
    knot2_path = [
        (0, 1, 2), (0, 2, 2), (1, 2, 2),
        (1, 2, 1), (1, 1, 1), (0, 1, 1),
        (0, 1, 2)
    ]

    # Knot 3
    knot3_path = [
        (2, 0, 1), (2, 0, 2), (2, 1, 2),
        (2, 1, 1), (3, 1, 1), (3, 0, 1),
        (2, 0, 1)
    ]
    
    # Calculate the length of each knot
    len1 = calculate_path_length(knot1_path)
    len2 = calculate_path_length(knot2_path)
    len3 = calculate_path_length(knot3_path)
    
    # Calculate the total length
    total_length = len1 + len2 + len3
    
    # Print the final equation
    print(f"The minimum total number of edges is the sum of the lengths of the three minimal components.")
    print(f"{len1} + {len2} + {len3} = {total_length}")

if __name__ == "__main__":
    main()