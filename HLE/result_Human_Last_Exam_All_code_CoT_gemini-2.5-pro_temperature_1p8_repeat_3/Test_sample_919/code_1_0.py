import sympy as sp

def solve_emi_shielding_force():
    """
    This function symbolically derives and prints the formula for the force per unit area
    on the conducting plane in the described EMI shielding application.
    """
    # Define symbolic variables
    mu_0, mu, K_0, a, y, d = sp.symbols('mu_0 mu K_0 a y d')
    
    # Based on the derivation, the y-component of the magnetic field H at x=d is:
    # H_y_at_d = (K_0 * sin(a*y)) / (cosh(a*d) + (mu_0/mu) * sinh(a*d))
    
    # The magnetic field B is mu_0 * H
    # B_y_at_d = mu_0 * H_y_at_d
    
    # The force per unit area (magnetic pressure) is given by B^2 / (2*mu_0), directed
    # along the negative x-axis.
    # Force_per_area = - (B_y_at_d**2) / (2*mu_0)
    # Force_per_area = - (mu_0**2 * H_y_at_d**2) / (2*mu_0)
    # Force_per_area = - (mu_0 / 2) * H_y_at_d**2
    
    # Numerator of the final expression
    numerator = - (mu_0 / 2) * K_0**2 * sp.sin(a*y)**2
    
    # Denominator of the final expression
    denominator = (sp.cosh(a*d) + (mu_0/mu) * sp.sinh(a*d))**2
    
    # The force is in the x-direction
    direction = "i_x"
    
    print("The formula for the force per unit area is derived as follows:")
    print("-" * 30)
    print(f"Force Direction: (-){direction}")
    # Using pretty print for mathematical expressions
    print("\nNumerator:")
    sp.pprint(numerator)
    print("\nDenominator:")
    sp.pprint(denominator)
    print("\nFinal Expression for Force per unit area (vector):")
    
    final_expression = numerator / denominator
    # For pretty printing the final fraction
    display_fraction = sp.Mul(numerator, 1/denominator, evaluate=False)

    # In Sympy, to display the vector, we can just print the scalar part
    # and mention the direction vector.
    sp.pprint(display_fraction)
    print(f"* ({direction})")
    
    print("\n" + "-"*30)
    print("This corresponds to option D.")


solve_emi_shielding_force()