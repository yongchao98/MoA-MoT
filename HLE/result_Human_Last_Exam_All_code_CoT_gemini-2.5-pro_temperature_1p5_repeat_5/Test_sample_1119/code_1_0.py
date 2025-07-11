import math

def solve_sequence_count():
    """
    Calculates the number of possible sequences based on N, K, and M.

    A sequence of K positive integers (a_1, ..., a_K) is valid if:
    1. 1 <= a_i <= N for all i.
    2. a_1 < a_2 < ... < a_K (strictly increasing).
    3. a_{i+1} - a_i <= M for all i.
    4. A given condition M*(K-1) < N holds.
    """
    
    # --- Parameters for the problem ---
    # N: The maximum value for any number in the sequence.
    N = 20
    # K: The length of the sequence.
    K = 5
    # M: The maximum increase between consecutive numbers.
    M = 4
    # ------------------------------------

    print(f"Finding the number of sequences for N={N}, K={K}, M={M}.")
    
    # The number of such sequences is given by the formula:
    # Sum_{s=0 to K-1} [(-1)^s * C(K-1, s) * C(N - s*M, K)]
    # where C(n, k) is the binomial coefficient "n choose k".

    def combinations(n, k):
        """Calculates combinations C(n, k), returns 0 if n < k or k < 0."""
        if k < 0 or n < k:
            return 0
        return math.comb(n, k)

    total_count = 0
    full_equation_str = []
    full_values_str = []

    for s in range(K):
        # Calculate combination values for the current term
        comb1_val = combinations(K - 1, s)
        n_for_comb2 = N - s * M
        comb2_val = combinations(n_for_comb2, K)

        # If the term's value is zero, skip it unless it's the first term.
        # This cleans up the output string.
        if comb1_val * comb2_val == 0 and s > 0:
            continue
        
        # Determine the sign for the inclusion-exclusion principle
        if s % 2 == 1: # Odd s -> subtract term
            term_val = -1 * comb1_val * comb2_val
            sign = " - "
        else: # Even s -> add term
            term_val = comb1_val * comb2_val
            sign = " + "
            
        total_count += term_val

        # Construct the string parts for this term to show the calculation
        term_equation_str = f"C({K - 1}, {s}) * C({n_for_comb2}, {K})"
        term_values_str = f"{comb1_val} * {comb2_val}"

        if s == 0:
            full_equation_str.append(term_equation_str)
            full_values_str.append(term_values_str)
        else:
            full_equation_str.append(sign + term_equation_str)
            full_values_str.append(sign + term_values_str)
            
    # Print the detailed calculation as requested
    print("\nCalculation based on the formula:")
    print(" ".join(full_equation_str))
    print("= " + " ".join(full_values_str))
    print(f"= {total_count}")
    
    # Final answer in the specified format is printed last
    print(f"\nFinal Answer: The number of possible sequences is {total_count}.")
    
    # Returning the value for the final tag
    return total_count

# Execute the function and capture the result for the final output tag
final_answer = solve_sequence_count()
print(f"<<<{final_answer}>>>")
