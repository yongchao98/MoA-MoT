import math

def combinations(n, k):
    """
    Calculates the binomial coefficient "n choose k".
    Returns 0 if k > n or k < 0.
    """
    if k < 0 or k > n:
        return 0
    # For large numbers, math.comb is preferred (Python 3.8+)
    # Using factorial for compatibility.
    if k == 0 or k == n:
        return 1
    if k > n // 2:
        k = n - k
    
    # This calculation is safe from intermediate overflow for reasonable inputs
    # and is generally faster than using math.factorial.
    res = 1
    for i in range(k):
        res = res * (n - i) // (i + 1)
    return res

def count_sequences(N, K, M):
    """
    Calculates the number of possible sequences based on the derived formula
    and prints the detailed calculation steps.

    Args:
        N: The maximum value for any number in the sequence.
        K: The length of the sequence.
        M: The maximum increase between consecutive numbers.
    """
    print(f"Calculating for N={N}, K={K}, M={M}")
    print("-" * 30)

    # Check the problem condition
    if not (M * (K - 1) < N):
        print("Warning: The condition M*(K-1) < N does not hold.")
        print(f"M*(K-1) = {M*(K-1)}, N = {N}")

    # Store parts of the equation strings to print them later
    eq_parts_symbolic = []
    eq_parts_evaluated_c1 = []
    eq_parts_fully_evaluated = []
    
    total = 0

    # The formula is a sum from i=0 to K-1
    for i in range(K):
        # Determine the sign of the term
        sign_val = (-1)**i
        
        # Calculate the first combination C(K-1, i)
        c1_n, c1_k = K - 1, i
        c1_val = combinations(c1_n, c1_k)
        
        # Calculate the second combination C(N - M*i, K)
        c2_n, c2_k = N - M * i, K
        c2_val = combinations(c2_n, c2_k)
        
        # Calculate the value of the current term in the sum
        term_val = sign_val * c1_val * c2_val
        total += term_val
        
        # --- Prepare the strings for printing the formula ---
        
        # Determine the sign string (+ or -) for nice printing
        if i == 0:
            sign_str = ""
        else:
            sign_str = " + " if sign_val > 0 else " - "
        
        # Remove the sign from the values for printing, as it's in sign_str
        abs_c1_val = c1_val
        abs_term_val = abs(c1_val * c2_val)

        eq_parts_symbolic.append(f"{sign_str}C({c1_n}, {c1_k}) * C({N} - {M}*{i}, {c2_k})")
        eq_parts_evaluated_c1.append(f"{sign_str}{abs_c1_val} * C({c2_n}, {c2_k})")
        eq_parts_fully_evaluated.append(f"{sign_str}{abs_c1_val} * {c2_val}")

    print("Formula: Sum_{i=0 to K-1} [ (-1)^i * C(K-1, i) * C(N - M*i, K) ]")
    print()
    print("Number of Sequences = " + "".join(eq_parts_symbolic))
    print("                    = " + "".join(eq_parts_evaluated_c1))
    print("                    = " + "".join(eq_parts_fully_evaluated))
    print(f"                    = {total}")
    return total

# --- Main execution ---
if __name__ == "__main__":
    # Example values for N, K, and M.
    # You can change these to solve for a different case.
    N = 20
    K = 5
    M = 4
    
    final_answer = count_sequences(N, K, M)
    # The final answer required by the problem format
    # print(f"\nFinal Answer: <<< {final_answer} >>>")