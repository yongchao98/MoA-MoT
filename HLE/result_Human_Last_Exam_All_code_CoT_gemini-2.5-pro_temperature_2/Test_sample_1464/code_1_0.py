import math

def find_and_sort_roots():
    """
    This function identifies the roots of the given polynomial based on algebraic analysis,
    calculates their numerical values, sorts them, and prints the result.
    """
    # From analyzing the polynomial's coefficients using Vieta's formulas,
    # the four roots are identified as sqrt(14), 2*sqrt(6), sqrt(34), and 2*sqrt(11).
    # We will now put them into a list with their symbolic names and numerical values.
    roots_data = [
        ("sqrt(14)", math.sqrt(14)),
        ("2*sqrt(6)", 2 * math.sqrt(6)),
        ("sqrt(34)", math.sqrt(34)),
        ("2*sqrt(11)", 2 * math.sqrt(11)),
    ]

    # Sort the roots based on their numerical value (the second element in each tuple).
    sorted_roots_data = sorted(roots_data, key=lambda item: item[1])

    # The prompt asks to "output each number in the final equation".
    # The final equation can be represented in factored form as (X - r1)(X - r2)(X - r3)(X - r4) = 0.
    # The "numbers" in this equation are the roots r1, r2, r3, r4.
    factored_parts = [f"(X - {r[0]})" for r in sorted_roots_data]
    factored_equation = " * ".join(factored_parts) + " = 0"
    
    print("The polynomial's final equation can be written in factored form:")
    print(factored_equation)

    print("\nThe four roots of the equation in increasing order are:")
    
    sorted_roots_numeric = []
    for symbolic_name, numeric_value in sorted_roots_data:
        sorted_roots_numeric.append(numeric_value)
        print(f"{symbolic_name:<10} (approximately {numeric_value:.10f})")

if __name__ == "__main__":
    find_and_sort_roots()