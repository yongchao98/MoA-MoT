import numpy as np

def solve_for_radius():
    """
    Calculates the radius of a circle tightly packing fourteen circles of radius one.
    The radius R is determined by finding the roots of a polynomial derived from the
    geometry of the optimal packing configuration.
    """
    
    # The polynomial is for x = R^2, where R is the radius of the large circle.
    # The coefficients of x^6 - 24x^5 + 204x^4 - 784x^3 + 1444x^2 - 1104x + 256 = 0
    coeffs = [1, -24, 204, -784, 1444, -1104, 256]

    print("The problem is to find the radius R of a circle that can pack 14 unit circles.")
    print("The radius R is a root of a known polynomial equation.")
    print("Let x = R^2. The equation for x is:")
    
    # Print the equation from the coefficients
    equation_str = ""
    for i, c in enumerate(coeffs):
        power = len(coeffs) - 1 - i
        # Format the term with its coefficient, variable, and power
        if c == 0:
            continue
        
        # Sign
        sign = " - " if c < 0 else " + "
        if i == 0:
            sign = "" if c > 0 else "-"
            
        # Coefficient
        coeff_val = abs(c)
        coeff_str = str(coeff_val) if coeff_val != 1 or power == 0 else ""

        # Variable and power
        var_str = f"x^{power}" if power > 1 else "x" if power == 1 else ""
        
        # For the constant term
        if power == 0:
            var_str = ""
        
        equation_str += f"{sign}{coeff_str}{var_str}"

    print(f"{equation_str} = 0\n")


    # Find the roots of the polynomial for x = R^2
    roots_x = np.roots(coeffs)
    
    # Filter for roots that are real and positive, as R^2 must be real and positive.
    possible_r_squared = [r.real for r in roots_x if np.isreal(r) and r.real > 0]
    
    # The largest root corresponds to the radius for this packing problem.
    if not possible_r_squared:
        print("No valid real solution found.")
        return
        
    R_squared = max(possible_r_squared)
    
    # Calculate R
    R = np.sqrt(R_squared)

    # Round the final answer to 4 significant digits.
    # We use a format specifier for this.
    R_rounded = float(f"{R:.4g}")

    print(f"The largest real root for R^2 is: {R_squared}")
    print(f"The corresponding radius R is: {R}")
    print("\n-------------------------------------------------")
    print(f"The radius of the circle up to 4 significant digits is: {R_rounded}")

solve_for_radius()