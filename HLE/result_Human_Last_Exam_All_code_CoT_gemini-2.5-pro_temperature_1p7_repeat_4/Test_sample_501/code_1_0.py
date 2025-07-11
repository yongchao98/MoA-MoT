import sympy

def solve_polymer_force():
    """
    This function derives and prints the force law for a thermally isolated
    freely jointed polymer chain.

    The derivation follows these steps:
    1.  The total entropy S = S_kin + S_conf is constant for a slow, adiabatic extension.
    2.  S_conf = C1 - (3 * k_B * x**2) / (2 * n * l**2)
    3.  S_kin = (N_dof / 2) * k_B * ln(K) + C2. For a large chain, N_dof is approx. 2n.
    4.  So, (n * k_B * ln(K)) - (3 * k_B * x**2) / (2 * n * l**2) = constant.
    5.  This implies K(x) = K(0) * exp((3 * x**2) / (2 * n**2 * l**2)).
        Let K(0) be the given initial kinetic energy, E(0).
    6.  The force F is the work done per unit extension, which for an adiabatic
        system is F = -dK/dx.
    7.  Differentiating K(x) with respect to x gives the force law.
    """
    # Define symbolic variables
    E0 = sympy.Symbol("E(0)") # Kinetic energy at zero extension
    x = sympy.Symbol("x")      # End-to-end separation
    n = sympy.Symbol("n")      # Number of mass points
    l = sympy.Symbol("ell")    # Length of each strut

    # Expression for Kinetic Energy K(x)
    # K(x) = E(0) * exp( (3 * x^2) / (2 * n^2 * l^2) )
    K = E0 * sympy.exp((3 * x**2) / (2 * n**2 * l**2))

    # The force is the negative derivative of the kinetic energy with respect to x
    F = -sympy.diff(K, x)

    # To ensure the numbers are explicitly shown, let's build the string manually.
    # The derived force is: F(x) = -(3 * E(0) * x) / (n^2 * l^2) * exp((3 * x^2) / (2 * n^2 * l^2))
    
    force_law_str = "F(x) = -(3 * E(0) * x / (n**2 * l**2)) * exp(3 * x**2 / (2 * n**2 * l**2))"
    
    print("The derived force law F(x) between the polymer ends is:")
    print(force_law_str)

solve_polymer_force()