import math

def count_sequences(N, K, M):
    """
    Calculates the number of possible sequences based on the given constraints.

    A sequence of K positive integers (a_1, ..., a_K) such that:
    1. 1 <= a_i <= N for all i
    2. a_1 < a_2 < ... < a_K
    3. a_{i+1} - a_i <= M for i=1,...,K-1
    4. M*(K-1) < N

    The number of sequences is calculated using the formula derived from
    inclusion-exclusion:
    Sum_{j=0}^{K-1} [ (-1)^j * C(K-1, j) * C(N - j*M, K) ]
    """

    # Check if K is positive, as the problem statement implies.
    if K <= 0:
        print("K must be a positive integer.")
        return

    # Check the given condition M(K-1) < N, which ensures N-jM >= 0 for j<=K-1
    if not (M * (K - 1) < N):
        print(f"The condition M(K-1) < N does not hold for N={N}, K={K}, M={M}.")
        print("The formula may not be applicable or may result in an error.")
        # Although the formula holds, the problem states this as a given.
        # We check it for correctness. math.comb would raise a ValueError if N-jM<0.
    
    print(f"For N={N}, K={K}, M={M}, the number of possible sequences is given by the formula:")
    print(f"Sum_{{j=0}}^{{{K-1}}} (-1)^j * C({K-1}, j) * C({N} - j*{M}, {K})\n")
    print("Calculation:")

    total_sequences = 0
    equation_symbolic_parts = []
    equation_numeric_parts = []
    
    # Summation from j=0 to K-1
    for j in range(K):
        sign_val = (-1)**j
        sign_str = "+ " if sign_val > 0 else "- "

        c1 = math.comb(K - 1, j)
        
        n_term = N - j * M
        # math.comb(n, k) returns 0 if k > n, which is what we want.
        # The condition M(K-1) < N ensures n_term >= 0.
        c2 = math.comb(n_term, K)
        
        term_value = sign_val * c1 * c2
        total_sequences += term_value

        # Build the string parts for printing the detailed equation
        equation_symbolic_parts.append(f"{sign_str}C({K-1},{j})*C({n_term},{K})")
        equation_numeric_parts.append(f"{sign_str}{c1}*{c2}")

    # Clean up the leading "+ " from the first term for display
    if equation_symbolic_parts[0].startswith("+ "):
        equation_symbolic_parts[0] = equation_symbolic_parts[0][2:]
    if equation_numeric_parts[0].startswith("+ "):
        equation_numeric_parts[0] = equation_numeric_parts[0][2:]

    # Print the detailed calculation steps
    print(" ".join(equation_symbolic_parts))
    print("= " + " ".join(equation_numeric_parts))
    print(f"= {total_sequences}")
    
    # This is for the final answer format as requested.
    return total_sequences

# --- Example Usage ---
# You can change these values to solve for a different case.
N_val = 20
K_val = 4
M_val = 5

final_answer = count_sequences(N_val, K_val, M_val)
# The final answer is wrapped in <<<>>> as requested.
# The code above has already printed the detailed breakdown.
# print(f"\n<<<{final_answer}>>>") 