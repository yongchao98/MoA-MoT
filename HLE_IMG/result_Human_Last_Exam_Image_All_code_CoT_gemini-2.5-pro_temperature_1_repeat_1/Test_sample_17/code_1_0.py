import sympy

def calculate_theoretical_stress_at_sharp_notch():
    """
    Calculates and explains the theoretical stress at the tip of a perfectly sharp wedge
    using symbolic mathematics, as per the principles of linear elasticity.
    """

    # Define the symbolic variables for the equation.
    # sigma_nom: The nominal uniaxial stress applied far from the notch.
    # d: The depth of the notch (a positive geometric constant).
    # rho: The radius of curvature at the tip of the notch.
    sigma_nom, d, rho = sympy.symbols('sigma_nom d rho', positive=True, real=True)

    # The formula for the stress concentration factor (Kt) for a notch is:
    # Kt = 1 + 2 * sqrt(d / rho)
    # The maximum stress (sigma_max) at the notch tip is Kt * sigma_nom.
    
    # We define the complete equation for the maximum stress.
    # Note: In the equation below, the numbers 1 and 2 are explicit.
    sigma_max_formula = (1 + 2 * sympy.sqrt(d / rho)) * sigma_nom

    print("The formula for the maximum theoretical stress (sigma_max) at the tip of a notch is:")
    print(f"sigma_max = (1 + 2 * sqrt(d / rho)) * sigma_nom")
    print("\nwhere:")
    print("  sigma_nom = Nominal applied stress")
    print("  d         = Depth of the notch")
    print("  rho       = Radius of curvature at the notch tip")

    print("\nFor a theoretically perfect sharp wedge, the tip radius 'rho' approaches 0.")
    print("We calculate the limit of the stress formula as rho -> 0 to find the stress at the tip.")
    
    # Calculate the limit of the stress formula as rho approaches 0 from the positive side.
    limit_stress = sympy.limit(sigma_max_formula, rho, 0, dir='+')

    print("\nResult:")
    print(f"limit(sigma_max) as rho->0 = {limit_stress}")
    print("\nConclusion: The theoretical stress at the tip of a perfectly sharp wedge is infinite.")

if __name__ == '__main__':
    calculate_theoretical_stress_at_sharp_notch()