import numpy as np

def solve_circle_packing_14():
    """
    Calculates the radius of a circle that tightly packs 14 unit circles.

    The solution is based on finding the roots of a known polynomial related
    to the optimal packing configuration for n=14. The variable 'x' in the
    polynomial represents the square of the radius (R^2) of the containing circle.

    The polynomial is:
    108*x^6 - 1224*x^5 + 4944*x^4 - 8984*x^3 + 8283*x^2 - 3618*x + 600 = 0
    """

    # Coefficients of the polynomial, from the constant term to the highest power
    coeffs = [600, -3618, 8283, -8984, 4944, -1224, 108]

    # Find all complex roots of the polynomial
    roots = np.polynomial.polynomial.polyroots(coeffs)

    # The square of the radius (R^2) must be a real, positive number.
    # We filter the roots to find all positive real roots.
    positive_real_roots = []
    for root in roots:
        if np.isreal(root) and root.real > 0:
            positive_real_roots.append(root.real)
    
    # For n=14, the radius of the minimal enclosing circle corresponds to the
    # largest of these positive real roots.
    if not positive_real_roots:
        print("Error: No positive real roots found.")
        return

    R_squared = max(positive_real_roots)

    # The radius R is the square root of R_squared.
    R = np.sqrt(R_squared)

    # The problem specifies that the radius of the small circles is 1.
    r = 1
    final_R = R * r

    print("For a packing of 14 circles, the squared radius (R^2) of the containing circle is the largest root of a specific polynomial.")
    print(f"The largest positive real root for R^2 is: {R_squared:.4f}")
    
    # Output the numbers used in the final equation per the instructions
    print(f"The radius is calculated as: sqrt({R_squared:.4f})")
    
    # Print the final result formatted to 4 significant digits
    print(f"The radius of the large circle is: {final_R:.4g}")

solve_circle_packing_14()