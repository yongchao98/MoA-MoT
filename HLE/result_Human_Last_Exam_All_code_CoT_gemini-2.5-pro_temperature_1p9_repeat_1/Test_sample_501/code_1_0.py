import numpy as np

def solve_polymer_force():
    """
    This script calculates and displays the force law for a thermally
    isolated freely jointed polymer chain.
    """
    # Define example parameters for the polymer model.
    n = 1000  # Number of mass points
    l = 1.0   # Length of one strut in nm
    x = 20.0  # Separation of the ends in nm
    E0 = 4.14e-21 # Kinetic energy at x=0 in Joules (approx. k_B * 300K)

    # The derived force law F(x) depends on the extension x,
    # the polymer properties n and l, and the initial kinetic energy E(0).
    # F(x) = - (2 * x * E(0) / (n^2 * l^2)) * exp(x^2 / (n^2 * l^2))

    # --- Calculation and Output ---

    # 1. Construct the string for the equation with numerical values
    equation_with_numbers = (
        f"F = - (2 * {E0:.2e} J * {x} m / ({n}^2 * {l}^2 m^2)) * "
        f"exp({x}^2 m^2 / ({n}^2 * {l}^2 m^2))"
    )
    
    print("The derived force law is:")
    print("F(x) = - (2 * E(0) * x / (n^2 * l^2)) * exp(x^2 / (n^2 * l^2))\n")
    
    print("For the example values:")
    print(f"  Number of segments (n): {n}")
    print(f"  Length of strut (l): {l} nm")
    print(f"  End-to-end separation (x): {x} nm")
    print(f"  Initial energy (E(0)): {E0:.2e} J\n")

    print("The equation with these values plugged in is:")
    # We use base SI units for the printout to be consistent, converting nm to m for E0 in J
    equation_with_si_units = (
        f"F = - (2 * {E0:.2e} * {x*1e-9:.1e} / ({n}^2 * ({l*1e-9:.1e})^2)) * "
        f"exp(({x*1e-9:.1e})^2 / ({n}^2 * ({l*1e-9:.1e})^2))"
    )
    print(equation_with_si_units)

    # 2. Calculate the final force value in Newtons
    n_squared_l_squared = n**2 * (l * 1e-9)**2 # in m^2
    x_m = x * 1e-9 # in m
    
    prefactor = - (2 * E0 * x_m) / n_squared_l_squared
    exponent = x_m**2 / n_squared_l_squared
    force_newtons = prefactor * np.exp(exponent)

    print(f"\nThe calculated restoring force is: {force_newtons:.3e} N")
    print("(A negative force indicates a restoring force, i.e., attraction between the ends)")

solve_polymer_force()