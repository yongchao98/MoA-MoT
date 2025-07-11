import math

def solve_sequence_count():
    """
    This function calculates the number of possible sequences based on the given constraints
    and prints the detailed calculation steps.
    """
    # --- Parameters ---
    # N: The maximum value for any number in the sequence.
    # K: The length of the sequence.
    # M: The maximum increase between consecutive numbers.
    N = 10
    K = 4
    M = 3
    # ------------------

    # The plan is to use the formula derived from the principle of inclusion-exclusion:
    # Result = Sum_{j=0}^{K-1} [ (-1)^j * C(K-1, j) * C(N - j*M, K) ]
    # where C(n, k) is the binomial coefficient "n choose k".
    
    print(f"Finding the number of sequences for N={N}, K={K}, M={M}:\n")
    
    total_count = 0
    
    # equation_parts will store the string components for each term in the sum
    # Format: [operator, symbolic_term, value_term, product_term]
    equation_parts = []
    
    # Loop through each term of the inclusion-exclusion formula
    for j in range(K):
        # Calculate the sign of the term: (-1)^j
        sign = (-1)**j
        
        # Calculate the first binomial coefficient: C(K-1, j)
        comb1_val = math.comb(K - 1, j)
        
        # Calculate the second binomial coefficient: C(N - j*M, K)
        n_comb2 = N - j * M
        k_comb2 = K
        
        # math.comb(n, k) returns 0 if n < k, which is the desired behavior.
        comb2_val = math.comb(n_comb2, k_comb2)

        # Calculate the full value of the term
        term_value = sign * comb1_val * comb2_val
        total_count += term_value
        
        # Determine the operator string ("+" or "-")
        op_str = " + " if sign > 0 else " - "
        # For the first term (j=0), omit the leading "+"
        if j == 0:
            op_str = "" if sign > 0 else "- "

        # Store the string representations of the current term
        equation_parts.append([
            op_str,
            f"C({K - 1}, {j})*C({N - j * M}, {K})",
            f"{comb1_val}*{comb2_val}",
            f"{abs(term_value)}"
        ])

    # Assemble and print the full equation from the stored parts
    eq_symbolic = "".join([p[0] + p[1] for p in equation_parts])
    eq_values = "".join([p[0] + p[2] for p in equation_parts])
    eq_products = "".join([p[0] + p[3] for p in equation_parts])
    
    print(f"Total = {eq_symbolic}")
    print(f"      = {eq_values}")
    print(f"      = {eq_products}")
    print(f"      = {total_count}")

# Execute the function
solve_sequence_count()