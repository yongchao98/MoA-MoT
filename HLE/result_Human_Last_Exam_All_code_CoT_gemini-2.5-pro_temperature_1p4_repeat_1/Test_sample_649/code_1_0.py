import sys

def print_lower_bound_formula():
    """
    This function prints the derived lower bound formula for the minimax risk R*_n,
    outputting each symbol or number of the final equation as requested.
    """
    # Using a list of strings to represent the formula components
    # to print them one by one as per the instruction.
    formula_parts = [
        "R*_n >= ", "1", "/", "2", " * ", "Φ", "(", "δ", " / ", "2", ")",
        " * ", "(", "1", " - ", "d_TV", "(", "P_0^n", ", ", "P", ")", ")"
    ]

    for part in formula_parts:
        # sys.stdout.write is used to avoid adding newlines automatically
        sys.stdout.write(part)
    
    # Add a final newline
    sys.stdout.write("\n\n")

    # Explanation of the terms in the formula
    explanation = [
        "where:",
        "  R*_n: The minimax risk for n samples.",
        "  Φ: The non-decreasing function from the loss definition.",
        "  δ: The minimum separation between the parameter θ(P_0) and the parameters θ(P_j) for j=1..N.",
        "  d_TV: The total variation distance between two probability distributions.",
        "  P_0^n: The product measure for n i.i.d. samples from distribution P_0.",
        "  P: The mixture distribution defined as P = (1/N) * sum_{j=1 to N} P_j^n."
    ]
    
    print("\n".join(explanation))


# Execute the function to print the result.
print_lower_bound_formula()
