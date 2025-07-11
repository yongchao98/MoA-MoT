import math

def solve_minimal_polynomial():
    """
    This function calculates and prints the minimal polynomial for the shortest 
    geodesic distance an ant can travel on a unit dodecahedron to return 
    to a vertex without passing through any others.

    The solution is based on the known squared distance L^2 = 13 + 6*sqrt(5).
    From this, we derive the minimal polynomial L^4 - 26*L^2 - 11 = 0.
    """

    # The coefficients of the minimal polynomial P(d) = a*d^4 + b*d^3 + c*d^2 + e*d + f = 0
    # The derived polynomial is d^4 - 26*d^2 - 11 = 0.
    coeffs = {
        4: 1,
        3: 0,
        2: -26,
        1: 0,
        0: -11
    }

    print("The minimal polynomial for the shortest distance 'd' is:")

    equation_parts = []
    # Iterate from the highest power down to 0 to build the equation string
    for power in sorted(coeffs.keys(), reverse=True):
        coeff = coeffs[power]
        
        # Skip terms with zero coefficient
        if coeff == 0:
            continue

        # Determine the sign
        if coeff > 0:
            sign = "+"
        else:
            sign = "-"
        
        # Get the absolute value of the coefficient
        coeff_abs = abs(coeff)

        # Format the variable part
        if power == 0:
            var_part = ""
        elif power == 1:
            var_part = "*d"
        else:
            var_part = f"*d^{power}"

        # Combine sign, coefficient, and variable part
        # e.g., "+ 1*d^4", "- 26*d^2", "- 11"
        if var_part:
            term = f"{sign} {coeff_abs}{var_part}"
        else:
            term = f"{sign} {coeff_abs}"
        
        equation_parts.append(term)
    
    # Join all parts and remove the leading "+ " if it exists
    final_equation = " ".join(equation_parts).lstrip("+ ")

    print(f"{final_equation} = 0")

solve_minimal_polynomial()
