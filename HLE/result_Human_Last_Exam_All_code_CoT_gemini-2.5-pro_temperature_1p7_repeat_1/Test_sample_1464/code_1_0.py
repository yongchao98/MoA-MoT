import sympy as sp

def find_polynomial_roots():
    """
    This function defines and solves the given polynomial equation symbolically
    and prints the roots in a structured format.
    """
    # Define the symbol X
    X = sp.Symbol('X')

    # Define the polynomial expression from the problem statement
    # The coefficients are grouped for clarity
    c3 = -(sp.sqrt(34) + sp.sqrt(14) + 2*sp.sqrt(11) + 2*sp.sqrt(6))
    c2 = (2*sp.sqrt(374) + 2*sp.sqrt(154) + 2*sp.sqrt(119) + 4*sp.sqrt(66) + 
          4*sp.sqrt(51) + 4*sp.sqrt(21))
    c1 = -(4*sp.sqrt(1309) + 4*sp.sqrt(714) + 8*sp.sqrt(561) + 8*sp.sqrt(231))
    c0 = 8*sp.sqrt(7854)

    poly = X**4 + c3*X**3 + c2*X**2 + c1*X + c0

    # Solve the equation poly = 0 for X
    try:
        roots = sp.solve(poly, X)
    except Exception as e:
        print(f"An error occurred during solving: {e}")
        return

    # Sort the roots in increasing order based on their numerical value
    sorted_roots = sorted(roots, key=lambda r: r.evalf())

    # To satisfy "output each number in the final equation!", we print the factored form
    factored_form = " * ".join([f"(X - ({sp.pretty(r)}))" for r in sorted_roots]) + " = 0"
    print("The factored form of the equation is:")
    print(factored_form)
    print("-" * 20)

    # Print the sorted roots clearly
    print("The 4 roots in increasing order are:")
    for i, root in enumerate(sorted_roots):
        print(f"Root {i+1}: {sp.pretty(root)}")

# Execute the function to find and print the roots
find_polynomial_roots()

<<<2*sqrt(11)>>>