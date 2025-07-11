import math

def solve_polynomial():
    """
    This function finds and prints the four roots of the given polynomial.
    """
    # Based on the analysis, the four roots are sqrt(14), 2*sqrt(6), sqrt(34), and 2*sqrt(11).
    # We will store them as a list of tuples, with their symbolic name and numeric value.
    roots = [
        ("sqrt(14)", math.sqrt(14)),
        ("2*sqrt(6)", 2 * math.sqrt(6)),
        ("sqrt(34)", math.sqrt(34)),
        ("2*sqrt(11)", 2 * math.sqrt(11))
    ]

    # Sort the list of roots based on their numeric value in increasing order.
    sorted_roots = sorted(roots, key=lambda x: x[1])

    print("The four roots of the polynomial in increasing order are:")
    # Print each root with its symbolic representation and numerical value.
    for i, (symbolic, numeric) in enumerate(sorted_roots):
        print(f"Root {i+1}: {symbolic} \u2248 {numeric:.4f}")

    # To satisfy the prompt "output each number in the final equation",
    # we display the polynomial in its factored form.
    print("\nThe final equation in factored form is:")
    equation_factors = [f"(X - {r[0]})" for r in sorted_roots]
    final_equation = " * ".join(equation_factors) + " = 0"
    print(final_equation)

# Execute the function
solve_polynomial()