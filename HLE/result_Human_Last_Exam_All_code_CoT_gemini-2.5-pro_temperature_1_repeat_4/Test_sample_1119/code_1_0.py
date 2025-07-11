import math

def combinations(n, k):
    """
    Calculates the binomial coefficient C(n, k), also known as "n choose k".
    Returns 0 if k > n or k < 0.
    """
    if k < 0 or k > n:
        return 0
    # math.comb is efficient and handles large numbers. Requires Python 3.8+
    return math.comb(n, k)

def solve_sequence_count():
    """
    Calculates the number of sequences of K positive integers up to N,
    where the sequence is strictly increasing and the gap between consecutive
    numbers is at most M.
    """
    # --- Parameters ---
    # You can change these values to solve for a different case.
    N = 20
    K = 5
    M = 4
    # ------------------

    print(f"Finding the number of sequences for N={N}, K={K}, M={M}")

    # Verification of the given condition M*(K-1) < N
    if not (M * (K - 1) < N):
        print(f"Error: The condition M(K-1) < N is not met.")
        print(f"Calculation: {M}*({K}-1) = {M*(K-1)}, which is not less than {N}.")
        return

    print("\nThe number of sequences is calculated using the formula:")
    print(f"  Sum_{{j=0 to {K-1}}} [ (-1)^j * C({K-1}, j) * C({N} - {M}*j, {K}) ]\n")

    print("--- Calculation Breakdown ---")
    
    total_sequences = 0
    final_eq_parts = []

    # The sum is from j=0 to K-1. The loop range(K) covers this.
    for j in range(K):
        sign = (-1)**j
        c1 = combinations(K - 1, j)
        n_for_c2 = N - M * j
        c2 = combinations(n_for_c2, K)
        
        term = sign * c1 * c2
        total_sequences += term

        # Build the string for the current term of the equation
        print(f"For j = {j}:")
        print(f"  Term = (-1)^{j} * C({K-1}, {j}) * C({N} - {M}*{j}, {K})")
        print(f"       = ({sign}) * C({K-1}, {j}) * C({n_for_c2}, {K})")
        print(f"       = ({sign}) * {c1} * {c2}")
        print(f"       = {term}\n")
        
        # Store term string for the final summary equation
        final_eq_parts.append(str(term))

    print("--- Final Result ---")
    # Creates a string like "15504 + -17472 + 4752 + -224 + 0" and cleans it up
    final_equation = " + ".join(final_eq_parts).replace("+ -", "- ")
    print(f"Total = {final_equation}")
    print(f"The total number of possible sequences is: {total_sequences}")
    
    return total_sequences

# Run the solver
final_answer = solve_sequence_count()
# The final answer will be printed in the <<<...>>> format below, as requested.
# The value is obtained from running the script.