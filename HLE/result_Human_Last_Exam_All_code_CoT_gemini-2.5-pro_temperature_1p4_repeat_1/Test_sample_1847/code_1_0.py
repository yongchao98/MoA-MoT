import numpy as np

def solve_shock_profile():
    """
    Calculates and prints the analytical solution for the density profile
    of a Mach sqrt(3) shock wave in a monatomic gas with Pr = 3/4.
    """

    # --- Step 1: Define physical constants ---
    # Prandtl number for the gas
    Pr = 3.0 / 4.0
    # Specific heat ratio for a monatomic ideal gas
    gamma = 5.0 / 3.0
    # Shock Mach number squared
    M2 = 3.0

    # --- Step 2: Calculate boundary conditions (Rankine-Hugoniot) ---
    # The density profile rho/rho_0 transitions from a normalized ambient
    # value of 1 to a post-shock value rho_1/rho_0.
    rho_0_norm = 1.0
    # The post-shock density ratio is given by the Rankine-Hugoniot jump condition:
    rho_1_norm = (gamma + 1) * M2 / (2 + (gamma - 1) * M2)

    # --- Step 3: Determine the coefficients for the tanh solution ---
    # The solution for the density profile has the form:
    # rho(x)/rho_0 = C1 + C2 * tanh((x/L)/D_prime)
    # C1 is the average density across the shock.
    C1 = (rho_0_norm + rho_1_norm) / 2.0
    # C2 is half the amplitude of the density change.
    C2 = (rho_1_norm - rho_0_norm) / 2.0

    # --- Step 4: Calculate the characteristic shock thickness ---
    # The denominator D_prime in the tanh argument determines the shock thickness.
    # It is derived from the transport properties (viscosity, thermal conductivity)
    # and gas properties. For this specific problem, it simplifies to:
    # D_prime = 4 * delta / L, where delta/L = (8 * Pr) / (3 * (gamma + 1))
    delta_over_L = (8.0 * Pr) / (3.0 * (gamma + 1.0))
    D_prime = 4.0 * delta_over_L

    # --- Step 5: Print the final analytical solution ---
    # The f-string formatting will output each number in the final equation.
    print("The analytical solution for the density profile is given by the equation:")
    print(f"rho(x)/rho_0 = {C1} + {C2} * tanh((x/L) / {int(round(D_prime))})")

solve_shock_profile()