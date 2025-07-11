import math

def solve_polynomial_roots():
    """
    This function finds and prints the four roots of the given polynomial in increasing order.
    The roots were identified by inspecting the polynomial's coefficients and applying Vieta's formulas.
    """
    
    # The four roots of the polynomial are sqrt(14), 2*sqrt(6), sqrt(34), and 2*sqrt(11).
    # We define them here as tuples containing their numerical value and symbolic representation.
    roots = [
        (math.sqrt(14), "sqrt(14)"),
        (2 * math.sqrt(6), "2*sqrt(6)"),
        (math.sqrt(34), "sqrt(34)"),
        (2 * math.sqrt(11), "2*sqrt(11)")
    ]

    # Sort the roots based on their numerical value in increasing order.
    roots.sort(key=lambda x: x[0])

    # Print the sorted roots.
    print("The 4 roots of the equation in increasing order are:")
    for i, (value, name) in enumerate(roots):
        print(f"Root {i+1}: {name} \u2248 {value:.7f}")

    # To fulfill the requirement "output each number in the final equation",
    # we print the equation in its factored form.
    factored_equation_parts = [f"(X - {name})" for value, name in roots]
    factored_equation = " * ".join(factored_equation_parts)
    print("\nThe equation in factored form is:")
    print(f"{factored_equation} = 0")


solve_polynomial_roots()

# <<<3.7416574, 4.8989795, 5.8309519, 6.6332496>>>