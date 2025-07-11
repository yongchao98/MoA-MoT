import sympy

def solve_dehn_twist_coefficient():
    """
    Calculates the fractional Dehn twist coefficient for (Da o Db)^9
    in the mapping class group of a torus with one boundary component.
    """

    # Step 1: Define matrices for Dehn twists Da and Db in SL(2, Z).
    # Da is a twist about the longitude 'a'.
    A = sympy.Matrix([[1, 1], [0, 1]])
    # Db is a twist about the meridian 'b'. We assume a standard basis where
    # the algebraic intersection number i(a, b) = 1.
    B = sympy.Matrix([[1, 0], [-1, 1]])

    # The matrix for the composition (Da o Db) is M = A * B.
    M = A * B

    # Step 2: Verify the boundary twist relation by finding the order of M.
    # We are looking for the smallest integer k such that M^k is the identity matrix.
    # This k corresponds to the power in the relation (Da o Db)^k = D_delta.
    k = 1
    M_k = M
    while not M_k.is_Identity:
        M_k = M_k * M
        k += 1

    # Step 3: Decompose the given map (Da o Db)^9.
    total_power = 9
    
    # The coefficient of the boundary twist D_delta is how many times (Da o Db)^k
    # fits into the total power.
    coefficient = total_power // k
    
    # The remainder is the power of the map that is left over.
    remainder_power = total_power % k

    # Step 4: Display the result.
    # The fractional Dehn twist coefficient is the integer 'coefficient'.
    print("The fundamental relation in the mapping class group of a punctured torus is:")
    print(f"(D_a \circ D_b)^{k} = D_\delta")
    print("\nUsing this relation, we can decompose the given map:")
    print(f"(D_a \circ D_b)^{total_power} = (D_a \circ D_b)^{k} \circ (D_a \circ D_b)^{remainder_power}")
    print(f"                 = (D_\delta)^{coefficient} \circ (D_a \circ D_b)^{remainder_power}")
    
    print("\nThe fractional Dehn twist coefficient is the exponent of D_delta in this decomposition.")
    print(f"The coefficient is: {coefficient}")
    
    print("\n--- Final Equation ---")
    print(f"(D_a \circ D_b)^{total_power} = (D_\delta)^{{{coefficient}}} \circ (D_a \circ D_b)^{{{remainder_power}}}")


solve_dehn_twist_coefficient()