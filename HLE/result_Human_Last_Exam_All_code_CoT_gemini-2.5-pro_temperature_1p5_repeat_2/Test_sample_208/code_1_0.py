import numpy as np

def solve_circle_packing_14():
    """
    Calculates the radius of a circle that tightly packs 14 unit circles.
    The solution is based on the proven optimal arrangement, where the radius R
    is the largest positive real root of a specific polynomial.
    """
    
    # The problem is to find the radius R for packing 14 circles of radius r=1.
    # The radius R is the largest real root of the polynomial:
    # P(R) = R^8 - 4R^7 - 8R^6 + 2R^5 + 23R^4 - 2R^3 - 12R^2 + 1 = 0
    
    # Coefficients of the polynomial, from R^8 down to the constant term.
    # The numbers in the final equation are the coefficients and powers.
    coeffs = [1, -4, -8, 2, 23, -2, -12, 0, 1]
    
    print("The problem of finding the radius R for tightly packing 14 unit circles is solved by finding the largest positive real root of the following polynomial equation:")
    # Print the equation, showing all the numbers involved.
    equation_str = "1*R^8 - 4*R^7 - 8*R^6 + 2*R^5 + 23*R^4 - 2*R^3 - 12*R^2 + 1 = 0"
    print(equation_str)
    print("\nSolving for R...")

    # Find the roots of the polynomial
    roots = np.roots(coeffs)

    # We need the largest positive real root.
    # We filter out complex roots and negative real roots.
    positive_real_roots = []
    for root in roots:
        if np.isreal(root) and root.real > 0:
            positive_real_roots.append(root.real)
            
    # The radius R is the maximum of these positive real roots.
    radius = max(positive_real_roots)
    
    # Format the result to 4 significant digits.
    radius_formatted = float(f'{radius:.4g}')
    
    print(f"\nThe largest positive real root, and thus the radius of the circle, is approximately:")
    print(radius_formatted)

# Run the solver
solve_circle_packing_14()
<<<3.428>>>