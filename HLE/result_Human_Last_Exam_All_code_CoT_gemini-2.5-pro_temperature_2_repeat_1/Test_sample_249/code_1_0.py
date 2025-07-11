import math

def solve():
    """
    Calculates the minimum possible diameter of a tree G with n+2 vertices and m leaves.
    """
    # Read n and m from user input
    try:
        n = int(input("Enter the value of positive integer n: "))
        m = int(input("Enter the value of positive integer m: "))
        if n <= 0 or m <= 0:
            print("n and m must be positive integers.")
            return
    except ValueError:
        print("Invalid input. Please enter integers.")
        return

    # A tree with at least 3 vertices must have at least 2 leaves.
    # Since n is a positive integer, n+2 >= 3, so m >= 2 must hold for a valid tree.
    if n + 2 >= 3 and m < 2:
        print(f"For a tree with {n+2} vertices, there must be at least 2 leaves, but m={m} was given.")
        return

    # I represents the number of internal vertices in the tree G.
    I = n + 2 - m

    # For a valid tree with positive n, the number of internal vertices must be at least 1.
    # m = n + 2 implies n = 0, which is not allowed.
    if I <= 0:
        print("For positive n, m must be less than n+2.")
        return
        
    min_diameter = 0
    if I == 1:
        # T_I is a single vertex. D(T_I) = 0. G is a star graph.
        min_diameter = 2
    elif I == 2:
        # T_I is a path of length 1. D(T_I) = 1.
        min_diameter = 3
    else:  # I >= 3
        # In this case, T_I must have at least 3 vertices.
        # A star graph on I vertices has I-1 leaves. This construction is
        # possible if the number of available leaves 'm' is at least 'I-1'.
        if m >= I - 1:
            # We can form T_I as a star, its diameter is 2.
            min_diameter = 4
        else:
            # Not enough leaves to form a star. We form a "spider" graph.
            # We distribute I-1 vertices into m legs around a central vertex.
            # This minimizes diameter for a fixed number of leaves for T_I.
            # We choose to have m leaves for T_I.
            I_minus_1 = I - 1
            q = I_minus_1 // m
            r = I_minus_1 % m
            
            # The diameter of T_I is the sum of the lengths of the two longest legs.
            if r == 0:
                # All m legs have length q. D(T_I) = q + q = 2q.
                min_diameter = 2 * q + 2
            elif r == 1:
                # One leg of length q+1, others of length q. D(T_I) = (q+1) + q = 2q+1.
                min_diameter = 2 * q + 1 + 2
            else:  # r >= 2
                # Two or more legs have length q+1. D(T_I) = (q+1) + (q+1) = 2q+2.
                min_diameter = 2 * q + 2 + 2

    print(f"For n = {n} and m = {m}:")
    print(f"The number of vertices is {n+2} and the number of leaves is {m}.")
    print(f"The number of internal vertices is |I| = n + 2 - m = {I}.")
    print(f"The minimum possible value for the diameter of G is {min_diameter}.")

solve()
