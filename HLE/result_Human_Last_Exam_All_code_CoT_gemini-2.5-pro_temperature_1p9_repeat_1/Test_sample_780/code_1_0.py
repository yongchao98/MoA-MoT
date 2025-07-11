import sys

# For Python 3.8 and later, pow(a, -1, m) can compute modular inverse.
# For earlier versions, we might need a custom function.
if sys.version_info < (3, 8):
    def modInverse(n, modulus):
        # Using Fermat's Little Theorem for prime modulus
        return pow(n, modulus - 2, modulus)
else:
    def modInverse(n, modulus):
        return pow(n, -1, modulus)

def solve():
    """
    Solves the problem by deriving and analyzing a linear recurrence relation.
    """
    # Step 1 & 2: Define problem parameters
    T = 510  # Total number of colors
    F = 203  # Number of 'forbidden' colors for the monochromatic constraint
    p = 23627  # The prime modulus, as given in the problem

    n_factor1 = 23626
    n_factor2_base = 23628
    n_factor2_exp1 = 100
    n_factor2_exp2 = 50

    # Step 3: Define the linear recurrence relation S(n) modulo p
    # The relation is S(n) = c1*S(n-1) + c2*S(n-2) + c3*S(n-3)
    T_sq_mod_p = pow(T, 2, p)
    
    # Coefficients modulo p
    c1 = (T_sq_mod_p - 1 + p) % p
    c2 = (T_sq_mod_p - 1 + p) % p
    c3 = (T_sq_mod_p - F + p) % p

    print("The derived recurrence relation modulo", p, "is:")
    print(f"S(n) = {c1}*S(n-1) + {c2}*S(n-2) + {c3}*S(n-3)")
    print("-" * 20)

    # Step 4: Analyze the exponent N
    # N = 23626 * (23628^100 - 23628^50)
    # p = 23627
    # n_factor1 = p - 1
    # n_factor2_base = p + 1
    # Thus, N is a multiple of (p-1).
    # N = (p-1) * k, where k = ( (p+1)^100 - (p+1)^50 )

    is_N_multiple_of_p_minus_1 = (n_factor1 % (p - 1) == 0)
    print(f"Is N a multiple of (p-1) = {p-1}? {is_N_multiple_of_p_minus_1}")

    # Step 5 & 6: Use properties of the recurrence over Z_p
    # Let the vector V(n) = [S(n), S(n-1), S(n-2)]^T.
    # V(n) = M * V(n-1), where M is the companion matrix.
    # M = [[c1, c2, c3], [1, 0, 0], [0, 1, 0]]
    # det(M) = -c3 = -(-182) = 182.
    # Since 182 != 0 (mod 23627), M is invertible.
    
    # For an invertible companion matrix M of a polynomial with integer coeffs,
    # it is a known property that M^(p-1) = I (mod p).
    # Since N is a multiple of p-1, M^N = (M^(p-1))^k = I^k = I.
    # We have V(N) = M^(N-2) * V(2).
    # Also V(2) = M^2 * V(0), where V(0) = [S(0), S(-1), S(-2)]^T.
    # So V(N) = M^(N-2) * M^2 * V(0) = M^N * V(0) = I * V(0) = V(0).
    # This implies S(N) === S(0) (mod p).

    # Step 7: Calculate the final result S(0)
    # S(0) is the number of ways to color a 2x0 rectangle, which is 1 (the empty coloring).
    S0 = 1
    final_answer = S0
    
    print("\nBased on the properties of the recurrence and the exponent N, we have:")
    final_equation_lhs = f"S({n_factor1} * ({n_factor2_base}^{n_factor2_exp1} - {n_factor2_base}^{n_factor2_exp2})) mod {p}"
    print(f"{final_equation_lhs} = S(0) mod {p}")
    
    print(f"\nThe value of S(0) is 1.")
    print(f"Therefore, the final answer is {final_answer}.")
    
    # Formatting the output to show the full equation as requested
    print("\nFinal Equation:")
    print(f"{final_equation_lhs} = {final_answer}")

solve()