import sympy

def solve_fdtc():
    """
    Calculates the fractional Dehn twist coefficient of (D_a o D_b)^9.
    """
    # Step 1: Define the matrices for the Dehn twists D_a and D_b.
    # On the first homology group H_1 of the torus, with basis {a, b}
    # such that the algebraic intersection i(a, b) = 1, the right-handed
    # Dehn twists D_a and D_b act as follows:
    # D_a: a -> a, b -> b + a
    # D_b: a -> a - b, b -> b
    # This gives the following matrices in SL(2,Z).
    M_a = sympy.Matrix([[1, 1], [0, 1]])
    M_b = sympy.Matrix([[1, 0], [-1, 1]])

    # Step 2: Compute the matrix M_phi for the mapping class phi = (D_a * D_b)^9.
    M_product = M_a * M_b
    # We calculate M_product^9. We can notice that M_product has order 6.
    # (M_product)^3 is -I, and (M_product)^6 is I.
    # So, (M_product)^9 = (M_product)^6 * (M_product)^3 = I * (-I) = -I.
    M_phi = M_product**9

    # Step 3: Determine N, the order of M_phi.
    # The matrix M_phi is -I. We need to find the smallest integer N > 0
    # such that (M_phi)^N = I.
    # (-I)^2 = I, so N = 2.
    N = 2

    # Step 4: Determine k.
    # The definition of the FDTC states that for our phi, phi^N = (D_delta)^k.
    # We have phi^N = phi^2 = ((D_a * D_b)^9)^2 = (D_a * D_b)^18.
    # We use the known relation in Mod(S_1,1): (D_a * D_b)^6 = (D_delta)^(-1).
    # So, (D_a * D_b)^18 = ((D_a * D_b)^6)^3 = ((D_delta)^(-1))^3 = (D_delta)^(-3).
    # This means k = -3.
    k = -3

    # Step 5: Calculate the final coefficient k/N.
    fdtc = sympy.Rational(k, N)

    print("The problem is to find the fractional Dehn twist coefficient of phi = (D_a o D_b)^9.")
    print("\nStep 1: The mapping class phi acts on homology via the matrix M_phi = (M_a * M_b)^9.")
    print(f"M_a = {M_a.tolist()}, M_b = {M_b.tolist()}")
    print(f"M_phi = {M_phi.tolist()}")
    
    print("\nStep 2: Find the order N of M_phi.")
    print(f"The matrix M_phi is -I, so its order is N = {N}.")

    print("\nStep 3: Relate phi^N to the boundary twist D_delta.")
    print(f"We have phi^N = phi^{N} = (D_a o D_b)^18.")
    print("Using the relation (D_a o D_b)^6 = (D_delta)^-1, we get:")
    print(f"phi^{N} = ((D_a o D_b)^6)^3 = ((D_delta)^-1)^3 = (D_delta)^-3.")
    print(f"So, we have the relation phi^N = (D_delta)^k with N = {N} and k = {k}.")

    print("\nStep 4: The fractional Dehn twist coefficient is k/N.")
    # The prompt requires outputting each number in the final equation.
    print(f"The final equation is: {k} / {N}")
    print(f"Result: {fdtc}")

solve_fdtc()