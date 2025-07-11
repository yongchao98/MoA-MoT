import sympy

def solve_factorization():
    """
    Calculates and verifies a left coprime factorization for the given transfer function H(s).
    """
    # Define the symbolic variable 's'
    s = sympy.Symbol('s')

    # Define the transfer function matrix H(s)
    H = sympy.Matrix([
        [(s - 1) / (s + 1), 1],
        [2 / (s**2 - 1), 0]
    ])

    # The derived left factorization matrices D(s) and N(s)
    D = sympy.Matrix([
        [1, s - 1],
        [0, s**2 - 1]
    ])

    N = sympy.Matrix([
        [1, 1],
        [2, 0]
    ])

    # --- Verification Step 1: Check if H = D^(-1) * N ---
    # Calculate D^(-1) * N
    D_inv_N = D.inv() * N
    # Simplify the resulting matrix to compare with H
    verification_H = sympy.simplify(D_inv_N)
    
    is_factorization_correct = (verification_H == H)

    # --- Verification Step 2: Check for left coprimeness ---
    # Form the augmented matrix M = [D, N]
    M = D.row_join(N)
    
    # The matrices are left coprime if the GCD of all 2x2 minors of M is a constant.
    minors = []
    # Generate all 2x2 minors
    for j1 in range(M.shape[1]):
        for j2 in range(j1 + 1, M.shape[1]):
            minor_matrix = M[:, [j1, j2]]
            minors.append(minor_matrix.det())
            
    # Calculate the GCD of all the minors
    if len(minors) > 0:
        current_gcd = minors[0]
        for i in range(1, len(minors)):
            current_gcd = sympy.gcd(current_gcd, minors[i])
    else:
        current_gcd = 0
        
    # Coprime if the GCD is a non-zero constant
    is_coprime = sympy.sympify(current_gcd).is_constant() and current_gcd != 0

    # --- Print the results ---
    print("The given transfer function H(s) is:")
    sympy.pprint(H)
    print("\nA left coprime factorization H(s) = D(s)^-1 * N(s) is given by:")
    
    print("\nD(s) =")
    sympy.pprint(D)
    
    print("\nN(s) =")
    sympy.pprint(N)
    
    print("\n--- Verification ---")
    print(f"1. Is D(s)^-1 * N(s) == H(s)?  {is_factorization_correct}")
    print(f"2. Are D(s) and N(s) left coprime? {is_coprime}")
    if is_coprime:
        print(f"   (The GCD of all 2x2 minors of [D(s) N(s)] is: {current_gcd})")

solve_factorization()

# The final answer format
final_D = "D(s) = [[1, s - 1], [0, s**2 - 1]]"
final_N = "N(s) = [[1, 1], [2, 0]]"
final_answer = f"{final_D}, {final_N}"
print(f"\n<<<{final_answer}>>>")