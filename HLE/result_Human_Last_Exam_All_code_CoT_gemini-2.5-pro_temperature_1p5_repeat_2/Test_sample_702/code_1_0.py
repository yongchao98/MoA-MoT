def solve_connectivity():
    """
    Calculates the connectivity of the map
    Sigma(Omega S^4 wedge Omega S^6) -> Omega(S^4 wedge S^6).
    """

    # The problem asks for the connectivity of the map
    # f: Sigma(Omega A wedge Omega B) -> Omega(A wedge B)
    # where A = S^4 and B = S^6.

    # Step 1: Define the dimensions of the spheres A and B.
    sphere_dim_A = 4
    sphere_dim_B = 6

    # Step 2: Determine the connectivity of spaces A and B.
    # The connectivity of a sphere S^n is (n-1).
    # This means pi_i(S^n) = 0 for all i < n.
    conn_A = sphere_dim_A - 1
    conn_B = sphere_dim_B - 1

    # Step 3: Apply the standard formula for the connectivity of the map.
    # The formula is conn(A) + conn(B) + 1.
    # This formula is valid because A and B are suspensions.
    final_connectivity = conn_A + conn_B + 1

    # Step 4: Output the result as a full equation, showing all numbers.
    print("The connectivity of the map is given by the formula conn(S^4) + conn(S^6) + 1.")
    print("conn(S^4) = 4 - 1 = 3")
    print("conn(S^6) = 6 - 1 = 5")
    print("So, the final calculation is:")
    print(f"{conn_A} + {conn_B} + 1 = {final_connectivity}")

solve_connectivity()