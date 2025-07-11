import math

def solve():
    """
    Calculates the minimum possible diameter of a tree with n+2 vertices and m leaves.
    """
    try:
        # In this environment, we can't get user input.
        # Using example values for n and m.
        # To use with other values, change them here.
        n_str = "10"
        m_str = "4"
        
        n = int(n_str)
        m = int(m_str)

        if n <= 0 or m <= 0:
            print("n and m must be positive integers.")
            return

        # A tree with n+2 vertices (where n>=1, so at least 3 vertices)
        # must have at least 2 leaves and at most n+1 leaves.
        if m < 2 or m > n + 1:
            print(f"For n={n}, the number of leaves m must be between 2 and {n+1}.")
            # For this problem, we assume the given n and m allow for a valid tree.
            # return

    except (ValueError, IndexError):
        print("Please provide two positive integers n and m.")
        return

    # Case 1: Vertex-centered construction ("spider")
    # Total vertices V = n + 2.
    # The sum of leg lengths is (V - 1) = n + 1.
    q = (n + 1) // m
    r = (n + 1) % m

    if r == 0:
        d_vertex = 2 * q
    elif r == 1:
        d_vertex = 2 * q + 1
    else: # r >= 2
        d_vertex = 2 * q + 2

    min_diameter = d_vertex
    d_edge = float('inf')

    # Case 2: Edge-centered construction ("double spider")
    # This is only possible if there are at least two internal nodes.
    # Number of internal nodes I = (n+2) - m. I >= 2 implies m <= n.
    if m <= n:
        min_edge_val = float('inf')
        # We partition m leaves into m1 and m2=m-m1.
        # We only need to check m1 from 1 to m//2 due to symmetry.
        for m1 in range(1, m // 2 + 1):
            m2 = m - m1
            
            # We need to partition n non-central vertices into V1 and V2=n-V1.
            # To minimize ceil(V1/m1) + ceil(V2/m2), we should choose V1
            # such that V1/m1 is close to V2/m2, which means V1 is close to n * m1 / m.
            # We only need to check the integers around this optimal value.
            v1_opt_float = n * m1 / m
            
            # Check a small range of integers around the fractional optimal value
            # This is an optimization over iterating all possible V1 values.
            # The range needs to be checked for validity: m1 <= V1 <= n-m2
            v1_candidates = {math.floor(v1_opt_float), math.ceil(v1_opt_float)}

            for v1 in v1_candidates:
                if m1 <= v1 <= n - m2:
                    v2 = n - v1
                    current_val = math.ceil(v1 / m1) + math.ceil(v2 / m2)
                    min_edge_val = min(min_edge_val, current_val)
        
        # If a valid partition was found, calculate the diameter
        if min_edge_val != float('inf'):
            d_edge = 1 + min_edge_val
            min_diameter = min(d_vertex, d_edge)

    print(f"For n = {n} and m = {m}:")
    print(f"The number of vertices is n+2 = {n+2}.")
    print(f"The number of leaves is m = {m}.")
    print(f"\nAnalyzing possible constructions to minimize diameter:")
    
    # Equation for vertex-centered diameter
    print(f"1. Vertex-centered diameter calculation:")
    print(f"   q = floor((n+1)/m) = floor(({n}+1)/{m}) = {q}")
    print(f"   r = (n+1) mod m = ({n}+1) mod {m} = {r}")
    if r == 0:
        print(f"   Diameter = 2*q = 2*{q} = {d_vertex}")
    elif r == 1:
        print(f"   Diameter = 2*q + 1 = 2*{q} + 1 = {d_vertex}")
    else:
        print(f"   Diameter = 2*q + 2 = 2*{q} + 2 = {d_vertex}")
    
    # Equation for edge-centered diameter
    if m <= n:
        print(f"2. Edge-centered diameter calculation:")
        print(f"   Diameter = min(1 + ceil(V1/m1) + ceil(V2/m2)) over partitions.")
        print(f"   Minimum found: {d_edge}")
    else:
        print(f"2. Edge-centered construction is not possible as m > n.")

    print(f"\nThe minimum possible value for the diameter is {min_diameter}.")

solve()