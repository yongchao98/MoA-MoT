import sympy as sp

def display_force_formula():
    """
    Constructs and prints the formula for the force per unit area
    on the conducting plane using sympy for clear representation.
    """
    # Define the symbolic variables
    mu_0 = sp.Symbol('mu_0')        # Permeability of free space
    mu = sp.Symbol('mu')            # Permeability of the magnetic material
    K_0 = sp.Symbol('K_0')          # Amplitude of the current sheet
    a = sp.Symbol('a')              # Spatial frequency of the current
    y = sp.Symbol('y')              # y-coordinate
    d = sp.Symbol('d')              # Thickness of the air gap

    # Numerator of the force expression
    # Term: (mu_0 / 2) * K_0^2 * sin(a*y)^2
    numerator = (mu_0 / 2) * K_0**2 * (sp.sin(a * y))**2

    # Denominator of the force expression
    # Term: [cosh(a*d) + (mu_0/mu)*sinh(a*d)]^2
    denominator = (sp.cosh(a * d) + (mu_0 / mu) * sp.sinh(a * d))**2

    # The full expression for the magnitude of the force per area
    force_magnitude = numerator / denominator

    # The vector expression for the force per unit area
    # The direction is along the x-axis (i_x)
    print("The force per unit y-z area on the x = d interface is:")
    
    # Print the equation string
    # We use sp.pretty_print for a more readable output
    ix_hat = sp.Symbol('i_x')
    force_vector_expr = force_magnitude * ix_hat
    
    # We build the string manually to match the desired format for "each number".
    # Since this is a formula, we present the components of the formula.
    print("\nvec{f}/area = (")
    print(f"  Numerator term 1 (Permeability of free space): {mu_0}")
    print(f"  Numerator term 2 (Constant): 1/2 = 0.5")
    print(f"  Numerator term 3 (Current Amplitude Squared): {K_0}**2")
    print(f"  Numerator term 4 (Spatial Variation): sin({a}*{y})**2")
    print(") / (")
    print(f"  Denominator term: [cosh({a}*{d}) + ({mu_0}/{mu})*sinh({a}*{d})]**2")
    print(") * i_x (unit vector in x direction)")


    print("\nFinal symbolic expression:")
    sp.pprint(force_vector_expr, use_unicode=True)

if __name__ == "__main__":
    display_force_formula()