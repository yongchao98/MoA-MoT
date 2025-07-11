import sympy as sp

def solve_emi_shielding_force():
    """
    This function prints the derived formula for the force per unit area on the conductor.
    The derivation involves solving Laplace's equation for the magnetic scalar potential
    and applying boundary conditions for a perfect conductor and a current sheet interface.
    The final force is calculated using the magnetic pressure formula.
    """
    # Define symbolic variables
    mu_0, mu, K_0, a, y, d = sp.symbols('mu_0 mu K_0 a y d', real=True, positive=True)
    x_hat = sp.Matrix([1, 0, 0])

    # Derived expression for the squared tangential magnetic field at x=d
    Denominator = sp.cosh(a * d) + (mu_0 / mu) * sp.sinh(a * d)
    Hy_squared = (K_0 * sp.sin(a * y) / Denominator)**2

    # Force per unit area vector
    # f_vector = - (mu_0 / 2) * Hy_squared * x_hat
    # Here, we print the components of the vector in a readable format.

    # Format the terms of the equation for clarity
    numerator = f"mu_0 * K_0**2 * sin(a*y)**2"
    denominator = f"[cosh(a*d) + (mu_0/mu) * sinh(a*d)]**2"
    direction = "i_x"

    print("The derived force per unit area on the conductor at x=d is:")
    print(f"f/area = - (1/2) * ({numerator}) / ({denominator}) * {direction}")

    # Explicitly show the components from the correct answer choice
    print("\nBreaking down the final equation from option D:")
    print("Term 1 (coefficient): -mu_0 / 2")
    print(f"Term 2 (numerator): K_0^2 * sin^2(ay)")
    print(f"Term 3 (denominator): [cosh(ad) + (mu_0/mu)*sinh(ad)]^2")
    print(f"Direction: x-hat (i_x)")
    
solve_emi_shielding_force()