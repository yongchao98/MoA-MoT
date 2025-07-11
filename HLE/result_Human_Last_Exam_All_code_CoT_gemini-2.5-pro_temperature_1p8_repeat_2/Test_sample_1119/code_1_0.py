import math

def solve_sequences():
    """
    Calculates the number of possible sequences based on the given constraints
    N, K, and M. It also prints the full calculation steps.
    """
    # --- Problem Parameters ---
    # N: The maximum value for any number in the sequence.
    # K: The length of the sequence.
    # M: The maximum increase between consecutive numbers.
    N = 10
    K = 4
    M = 3
    # --------------------------

    print(f"Calculating for N={N}, K={K}, M={M}")
    print(f"The given condition M*(K-1) < N is: {M*(K-1)} < {N}, which is {M*(K-1) < N}.")
    if not M*(K-1) < N:
        print("Warning: The condition M*(K-1) < N is not met.")

    def combinations(n, k):
        """Helper function for combinations, returns 0 if n < k."""
        if k < 0 or n < k:
            return 0
        return math.comb(n, k)

    total_sequences = 0
    calculation_parts = []

    # Loop through the sum based on the inclusion-exclusion formula
    for s in range(K):
        # Calculate combination terms
        c1 = combinations(K - 1, s)
        arg_c2 = N - s * M
        c2 = combinations(arg_c2, K)
        
        term_value = c1 * c2

        if s % 2 == 1:
            # Odd 's' corresponds to a subtraction term
            sign = -1
        else:
            # Even 's' corresponds to an addition term
            sign = 1
        
        total_sequences += sign * term_value
        
        # Store parts for detailed printing
        calculation_parts.append({
            "sign_str": " - " if sign == -1 else " + ",
            "eq_str": f"C({K-1}, {s}) * C({arg_c2}, {K})",
            "val_str": f"{c1} * {c2}",
            "term_val_str": str(term_value)
        })

    # Print the full equation and calculation steps
    print("\nThe number of sequences is given by the formula:")
    
    # Build equation string
    eq_line = []
    vals_line = []
    calc_line = []
    
    for i, part in enumerate(calculation_parts):
        sign = part["sign_str"]
        if i == 0 and sign == " + ":
            sign = ""
        
        eq_line.append(f'{sign}{part["eq_str"]}')
        vals_line.append(f'{sign}{part["val_str"]}')
        calc_line.append(f'{sign}{part["term_val_str"]}')
        
    print(" ".join(eq_line))
    print("= " + " ".join(vals_line))
    print("= " + " ".join(calc_line))
    print(f"= {total_sequences}")
    
    # Return the final computed answer to be captured
    return total_sequences

# Run the function and capture the answer
final_answer = solve_sequences()
# The required final output format
# print(f"<<<{final_answer}>>>")