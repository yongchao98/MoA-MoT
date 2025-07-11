import numpy as np

def calculate_nc_mass_change():
    """
    Calculates the percentage change in the effective mass of a Dirac nucleon
    due to noncommutative quantum mechanics.
    """
    # Step 1: Define constants and convert to consistent units (MeV)
    M = 938.93  # Nucleon mass (MeV)
    Cs2 = 267.1  # Scalar coupling constant (dimensionless)
    kF_inv_fm = 1.42  # Fermi momentum (fm^-1)
    R_fm = 5.5  # Radius of heavy nuclei (fm)
    M_star_ratio = 0.78  # Effective mass ratio M*/M
    nu = 4  # Degeneracy factor for nuclear matter

    # Conversion factor in natural units (ħ=c=1)
    hbar_c = 197.3  # MeV fm

    # Derived quantities in MeV units
    M_star = M_star_ratio * M
    kF = kF_inv_fm * hbar_c
    R = R_fm / hbar_c # R is in MeV^-1

    print("--- Input Parameters (in MeV units) ---")
    print(f"Nucleon mass M = {M:.2f} MeV")
    print(f"Effective mass M* = {M_star:.2f} MeV")
    print(f"Scalar coupling constant C_s^2 = {Cs2}")
    print(f"Fermi momentum k_F = {kF:.2f} MeV")
    print(f"Nuclear radius R = {R:.5f} MeV^-1")
    print("-" * 20)

    # Step 2 & 3: Calculate the value of the momentum integral part
    # The integral to solve is J = Integral from 0 to kF of [k^3 / (k^2 + M*^2)^(3/2)] dk
    # The solved definite integral is:
    # J_val = [sqrt(k^2 + M*^2) + M*^2 / sqrt(k^2 + M*^2)] evaluated from 0 to kF
    sqrt_term_kF = np.sqrt(kF**2 + M_star**2)
    val_at_kF = sqrt_term_kF + (M_star**2 / sqrt_term_kF)
    val_at_0 = np.sqrt(0**2 + M_star**2) + (M_star**2 / np.sqrt(0**2 + M_star**2))
    J = val_at_kF - val_at_0
    
    print("--- Intermediate Calculations ---")
    print(f"Value of the momentum integral J = {J:.4f} MeV")

    # The full change in mass is ΔM* = η * F, where F is the rest of the expression.
    # We first calculate F = ΔM*/η.
    # From simplifying the formula, we have:
    # F = (Cs2 / M^2) * M* * J * (3 * nu * R) / (32 * pi^2)
    delta_M_star_over_eta = (Cs2 / M**2) * M_star * J * R * (3 * nu / (32 * np.pi**2))
    print(f"The change in mass per unit η (ΔM*/η) = {delta_M_star_over_eta:.6f} MeV^-1")

    # Step 4: Assume η = 1/R^2 based on the characteristic length scale of the system.
    eta = (1/R)**2
    print(f"Assuming η = 1/R^2 = {eta:.2f} MeV^2")

    # Step 5: Calculate the final mass change and percentage change.
    delta_M_star = delta_M_star_over_eta * eta

    percentage_change = (delta_M_star / M_star) * 100
    
    print("-" * 20)
    print("--- Final Calculation ---")
    print(f"The total change in effective mass is ΔM* = {delta_M_star:.4f} MeV")
    print(f"The percentage change is (ΔM* / M*) * 100%")
    print(f"Percentage Change = ({delta_M_star:.4f} MeV / {M_star:.2f} MeV) * 100% = {percentage_change:.4f}%")
    print(f"This is a positive change, approximately +{percentage_change:.2f}%.")


calculate_nc_mass_change()