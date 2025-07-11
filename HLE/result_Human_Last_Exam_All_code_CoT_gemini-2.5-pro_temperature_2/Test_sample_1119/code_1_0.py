import math

def calculate_sequences(N, K, M):
    """
    Calculates the number of possible sequences based on the derived formula.

    The number of sequences is given by the formula:
    Sum_{s=0 to K-1} [(-1)^s * C(K-1, s) * C(N - s*M, K)]
    where C(n, k) is the binomial coefficient "n choose k".
    """
    # Check if the problem's condition is met. This is for context.
    if M * (K - 1) >= N:
        print(f"Warning: The condition M*(K-1) < N is not met for N={N}, K={K}, M={M}.")
        print("The formula may not apply or the problem definition might be ill-posed.")
        # Although the formula might still yield a result, the context from the prompt is violated.
        # We'll proceed with calculation but this warning is important.

    total_count = 0
    formula_str_parts = []
    values_str_parts = []
    
    # The sum goes from s = 0 to K-1.
    for s in range(K):
        sign = (-1)**s
        
        # First binomial coefficient: C(K-1, s)
        try:
            comb1 = math.comb(K - 1, s)
        except ValueError:
            # This happens if s > K-1, loop condition should prevent this.
            comb1 = 0

        # Second binomial coefficient: C(N - s*M, K)
        n_for_comb2 = N - s * M
        comb2 = 0
        # math.comb(n,k) requires n >= k. If n < k, the combination is 0.
        if n_for_comb2 >= K:
            comb2 = math.comb(n_for_comb2, K)
            
        # Calculate the term for this value of s
        term = sign * comb1 * comb2
        total_count += term
        
        # Build the strings for printing the formula and values
        sign_str_op = " - " if sign == -1 else " + "
        # First term doesn't need a leading sign
        if s == 0:
            sign_str_op = ""
        
        formula_str_parts.append(f"{sign_str_op}C({K-1}, {s})*C({n_for_comb2}, {K})")
        values_str_parts.append(f"{sign_str_op}{comb1}*{comb2}")
        
    print(f"For N={N}, K={K}, M={M}, the number of sequences is calculated as:")
    # Print the formula with values
    print("".join(formula_str_parts))
    # Print the evaluated terms
    print("= " + "".join(values_str_parts).lstrip(" +"))
    # Print the final sum
    print(f"= {total_count}")
    return total_count

# --- Main execution ---
# Example values (you can change these)
N = 10
K = 4
M = 3

# Calculate and print the result
final_answer = calculate_sequences(N, K, M)