from sage.all import EllipticCurve, DirichletGroup, ComplexField

def solve_elliptic_curve_problem():
    """
    Solves the problem by calculating the rank and L-function values
    for the given elliptic curve and characters.
    """
    # Step 1: Define the Elliptic Curve
    # The equation is y^2 + y = x^3 - x^2 - 10x - 20.
    # This corresponds to coefficients [a1, a2, a3, a4, a6] = [0, -1, 1, -10, -20].
    try:
        E = EllipticCurve([0, -1, 1, -10, -20])
    except Exception as e:
        print(f"Error creating elliptic curve: {e}")
        return

    # Step 2: Compute the rank r
    r = E.rank()

    # Step 3: Define the Dirichlet Characters
    # We need the two primitive cubic characters of conductor 7.
    # We use a complex field with sufficient precision for numerical calculations.
    try:
        G = DirichletGroup(7, ComplexField(100))
        cubic_chars = [chi for chi in G if chi.order() == 3]
        if len(cubic_chars) < 2:
            print("Could not find two cubic characters of conductor 7.")
            return
        chi1 = cubic_chars[0]
        chi2 = cubic_chars[1]
    except Exception as e:
        print(f"Error creating Dirichlet characters: {e}")
        return

    # Step 4: Compute the leading coefficients a and b
    # The central_value() method computes L^(k)(E, 1, chi)/k! where k is the
    # order of vanishing, which is the leading Taylor coefficient.
    try:
        a = E.lseries().twist(chi1).central_value()
        b = E.lseries().twist(chi2).central_value()
    except Exception as e:
        print(f"Error computing L-function central values: {e}")
        return

    # Step 5: Calculate the final sum and print the results
    total = r + a + b

    # Print the equation with each number as requested
    print(f"The rank r is: {r}")
    print(f"The leading coefficient a is: {a}")
    print(f"The leading coefficient b is: {b}")
    print("\nThe final equation is r + a + b:")
    # The sum a+b is real, so we can represent the total cleanly.
    print(f"{r} + ({a}) + ({b}) = {total}")

    # Round the final result to four decimal places
    final_answer = round(float(total.real), 4)
    print(f"\nThe value of r + a + b rounded to four decimal places is: {final_answer}")

# Execute the function
solve_elliptic_curve_problem()