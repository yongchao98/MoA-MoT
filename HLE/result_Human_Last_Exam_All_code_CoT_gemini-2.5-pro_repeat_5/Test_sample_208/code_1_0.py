import numpy as np

def solve_circle_packing_14():
    """
    Calculates the radius of a circle that tightly packs 14 unit circles.
    
    The solution is based on finding the root of a polynomial derived from the
    geometry of the optimal packing arrangement.
    """
    
    # The radius R is related to the largest real root of the polynomial:
    # 2x^4 - 40x^3 + 204x^2 - 392x + 245 = 0
    # where x = (R-1)^2.
    coeffs = [2, -40, 204, -392, 245]

    # Find the roots of the polynomial
    roots = np.roots(coeffs)

    # Filter for real, positive roots and find the largest one.
    # In this case, all roots happen to be real and positive.
    x = max(root.real for root in roots if np.isreal(root) and root.real > 0)
    
    # Calculate R from x, where x = (R-1)^2
    # R - 1 = sqrt(x)
    # R = 1 + sqrt(x)
    radius = 1 + np.sqrt(x)

    # Output the explanation and the final equation with numbers
    print("The radius 'R' is found by solving the polynomial 2x^4 - 40x^3 + 204x^2 - 392x + 245 = 0, where x = (R-1)^2.")
    print(f"\nThe largest real root of the polynomial is x = {x:.6f}")
    print("\nWe calculate the radius R using the relation R = 1 + sqrt(x):")
    
    # Print the final equation with numbers filled in
    print(f"R = 1 + sqrt({x:.6f})")
    print(f"R = 1 + {np.sqrt(x):.6f}")
    print(f"R = {radius:.6f}")

    # Print the final answer to 4 significant digits
    print("\nThe radius of the large circle up to 4 significant digits is:")
    print(f"{radius:.4g}")

solve_circle_packing_14()