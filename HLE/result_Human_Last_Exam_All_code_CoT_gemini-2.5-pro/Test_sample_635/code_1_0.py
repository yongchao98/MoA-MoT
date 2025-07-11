import sympy

def calculate_left_coprime_factorization():
    """
    Calculates and verifies a left coprime factorization for the given
    transfer function matrix H(s).
    """
    # Define the symbolic variable 's'
    s = sympy.Symbol('s')

    # Define the transfer function matrix H(s)
    H_s = sympy.Matrix([
        [(s - 1) / (s + 1), 1],
        [2 / (s**2 - 1), 0]
    ])
    print("Given transfer function H(s):")
    sympy.pprint(H_s)
    print("-" * 30)

    # Step 1 & 2: Propose a polynomial matrix D(s).
    # Based on the condition D(s)H(s) = N(s), where N(s) must be a polynomial matrix,
    # we can derive the structure of D(s).
    # A minimal D(s) can be found by making it upper triangular and satisfying the polynomial conditions.
    # d21=0. For d22*h21 to be polynomial, d22 must be a multiple of s**2-1. Simplest choice: d22 = s**2-1.
    # For d11*h11 + d12*h21 to be polynomial, d11(s-1)^2+2*d12 must be divisible by s**2-1.
    # Simplest choice: d11=1, d12=s-1.
    
    D_s = sympy.Matrix([
        [1, s - 1],
        [0, s**2 - 1]
    ])
    
    print("Proposed polynomial matrix D(s):")
    sympy.pprint(D_s)
    print("This matrix is chosen to be as simple as possible while ensuring D(s)H(s) results in a polynomial matrix N(s).")
    print("-" * 30)

    # Step 3: Calculate N(s) = D(s) * H(s)
    # Sympy's simplify() is used to resolve the expressions into their simplest polynomial form.
    N_s = sympy.simplify(D_s * H_s)

    print("Calculated polynomial matrix N(s) = D(s) * H(s):")
    sympy.pprint(N_s)
    print("-" * 30)
    
    # Step 4: Verification
    print("Verification steps:")
    # a) Verify that H(s) = D(s)^-1 * N(s)
    try:
        H_check = sympy.simplify(D_s.inv() * N_s)
        is_correct_factorization = (H_check == H_s)
        print(f"1. Verification of H(s) = D(s)^-1 * N(s): {is_correct_factorization}")
    except Exception as e:
        print(f"1. Verification failed: {e}")

    # b) Verify that D(s) and N(s) are left coprime
    # This is true if the rank of [D(s) N(s)] is 2 for all complex s.
    # This is equivalent to the greatest common divisor (GCD) of all 2x2 minors of [D(s) N(s)] being a constant.
    M_s = D_s.row_join(N_s)
    print("\n2. Verifying coprimeness of D(s) and N(s).")
    print("   The concatenated matrix [D(s) N(s)] is:")
    sympy.pprint(M_s)
    
    minors = []
    # There are 6 possible 2x2 minors in a 2x4 matrix.
    for j1 in range(M_s.shape[1]):
        for j2 in range(j1 + 1, M_s.shape[1]):
            minor_matrix = M_s[:, [j1, j2]]
            minors.append(sympy.simplify(minor_matrix.det()))
    
    # Remove zero minors for clarity
    non_zero_minors = [m for m in minors if m != 0]
    print("\n   The 2x2 minors of [D(s) N(s)] are:")
    sympy.pprint(non_zero_minors)
    
    # Calculate the GCD of all the minors
    if len(non_zero_minors) > 0:
        current_gcd = non_zero_minors[0]
        for i in range(1, len(non_zero_minors)):
            current_gcd = sympy.gcd(current_gcd, non_zero_minors[i])
    else: # Should not happen if rank is > 0
        current_gcd = 0
        
    print(f"\n   The greatest common divisor (GCD) of the minors is: {current_gcd}")

    is_coprime = sympy.sympify(current_gcd).is_constant() and current_gcd != 0
    print(f"   Since the GCD is a non-zero constant, the rank of [D(s) N(s)] is 2 for all s.")
    print(f"   Therefore, D(s) and N(s) are left coprime: {is_coprime}")
    print("-" * 30)
    
    print("Final result for the left coprime factorization H(s) = D(s)^-1 * N(s):")
    print("\nThe matrix D(s) is:")
    print(f"Row 1: [ {D_s[0,0]}, {D_s[0,1]} ]")
    print(f"Row 2: [ {D_s[1,0]}, {D_s[1,1]} ]")
    
    print("\nThe matrix N(s) is:")
    print(f"Row 1: [ {N_s[0,0]}, {N_s[0,1]} ]")
    print(f"Row 2: [ {N_s[1,0]}, {N_s[1,1]} ]")

calculate_left_coprime_factorization()