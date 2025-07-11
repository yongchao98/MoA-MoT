import numpy as np

def calculate_mass_change():
    """
    Calculates the percentage change in the effective mass of a Dirac nucleon
    due to non-commutative geometry.
    """
    # Constants provided in the problem
    M = 938.93  # MeV (Nucleon mass)
    Cs2 = 267.1  # Dimensionless (Scalar coupling constant)
    kF_inv_fm = 1.42  # fm^-1 (Fermi momentum)
    R_fm = 5.5  # fm (Radius of heavy nuclei)
    nu = 4  # Degeneracy factor
    M_star_ratio = 0.78  # M*/M

    # Conversion factor in natural units (hbar = c = 1)
    hc = 197.327  # MeV fm

    # Convert all quantities to a consistent unit system (MeV)
    kF = kF_inv_fm * hc  # MeV
    R = R_fm / hc  # MeV^-1
    M_star = M_star_ratio * M  # MeV

    # --- Step 1: Calculate the term in the square brackets ---
    kF2 = kF**2
    M_star2 = M_star**2
    sqrt_term = np.sqrt(kF2 + M_star2)
    bracket_term = (kF2 + 2 * M_star2) / sqrt_term - 2 * M_star

    # --- Step 2: Calculate the coefficient of eta in the relative mass change ---
    # The formula for relative change is dM*/M* = eta * (coefficient)
    # dM_star / eta = (Cs2 / M**2) * (3 * nu * M_star * R) / (32 * np.pi**2) * bracket_term
    # (dM_star / M_star) / eta = (Cs2 / M**2) * (3 * nu * R) / (32 * np.pi**2) * bracket_term
    
    delta_M_star_over_eta = (Cs2 / M**2) * (3 * nu * M_star * R) / (32 * np.pi**2) * bracket_term
    
    rel_change_per_eta = delta_M_star_over_eta / M_star

    # --- Step 3: Estimate eta and calculate the final percentage change ---
    # Assume eta ~ 1/R^2 based on physical scaling arguments (eta*theta ~ 1, theta ~ R^2)
    eta = 1 / R**2
    
    # Calculate the total relative change
    relative_change = rel_change_per_eta * eta
    
    # Convert to percentage
    percentage_change = relative_change * 100

    # --- Step 4: Print the results ---
    print("--- Input Parameters ---")
    print(f"Nucleon Mass (M): {M:.2f} MeV")
    print(f"Effective Mass (M*): {M_star:.2f} MeV")
    print(f"Fermi Momentum (k_F): {kF:.2f} MeV")
    print(f"Nuclear Radius (R): {R_fm} fm = {R:.4f} MeV^-1")
    print(f"Scalar Coupling (Cs^2): {Cs2}")
    
    print("\n--- Intermediate Calculations ---")
    print(f"Value of the bracket term: {bracket_term:.4f} MeV")
    print(f"Relative change per unit of eta: {rel_change_per_eta:.4e} MeV^-2")
    print(f"Estimated value for eta = 1/R^2: {eta:.2f} MeV^2")

    print("\n--- Final Calculation ---")
    print("Percentage Change = (Relative Change per unit eta) * eta * 100%")
    print(f"Final Equation: {percentage_change:.3f}% = {rel_change_per_eta:.4e} * {eta:.2f} * 100")
    
    print(f"\nCalculated Percentage Change: +{percentage_change:.3f}%")

calculate_mass_change()
<<<A>>>