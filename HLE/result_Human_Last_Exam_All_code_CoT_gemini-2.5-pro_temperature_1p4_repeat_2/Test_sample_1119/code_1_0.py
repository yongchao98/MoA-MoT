import math

def solve_and_print_sequence_count(N, K, M):
    """
    Calculates the number of possible sequences based on the problem description
    and prints the detailed calculation.

    A sequence a_1, ..., a_K of positive integers satisfies:
    1. 1 <= a_i <= N
    2. a_1 < a_2 < ... < a_K
    3. a_{i+1} - a_i <= M
    4. M*(K-1) < N

    The number of such sequences can be found using the formula:
    Sum_{j=0}^{K-1} (-1)^j * C(K-1, j) * C(N - M*j, K)
    """

    # Helper for combinations, returns 0 if n < k
    def combinations(n, k):
        if k < 0 or n < k:
            return 0
        return math.comb(n, k)

    # --- Print Header ---
    print("Problem Parameters:")
    print(f"N (max value in sequence): {N}")
    print(f"K (length of sequence): {K}")
    print(f"M (max increase between elements): {M}")
    print("-" * 30)
    
    # --- Check Condition ---
    if M * (K - 1) >= N:
        print(f"Warning: The condition M*(K-1) < N is not met for the given values.")
        print("The formula is derived based on this condition holding.")
        # The calculation will proceed, but the result might not be meaningful
        # for the original problem context if this condition is violated.
    
    # --- Calculation ---
    symbolic_parts = []
    value_parts = []
    total_sequences = 0

    for j in range(K):
        sign = (-1)**j
        sign_char = "+" if sign > 0 else "-"

        # Symbolic Part
        symbolic_term = f"C({K-1}, {j})*C({N-M*j}, {K})"
        if j > 0:
            symbolic_parts.append(f" {sign_char} ")
        symbolic_parts.append(symbolic_term)

        # Value Part & Total Calculation
        comb1 = combinations(K-1, j)
        comb2 = combinations(N-M*j, K)
        
        value_term = f"{comb1}*{comb2}"
        if j > 0:
            value_parts.append(f" {sign_char} ")
        value_parts.append(value_term)
        
        total_sequences += sign * comb1 * comb2

    # --- Print Final Equation ---
    final_equation_string = "".join(symbolic_parts) + " = " + "".join(value_parts) + f" = {total_sequences}"
    
    print("The calculation based on the formula is:")
    print(final_equation_string)

# --- Main execution ---
if __name__ == "__main__":
    # Example values for N, K, and M.
    # You can change these to solve for different parameters.
    N = 20
    K = 5
    M = 4
    
    solve_and_print_sequence_count(N, K, M)