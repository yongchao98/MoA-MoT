import math

def calculate_mass_change():
    """
    Calculates the percentage change in the effective mass of the Dirac nucleon 
    due to non-commutative geometry.
    """
    # --- Given constants ---
    # Nucleon mass in MeV
    M = 938.93
    # Scalar coupling constant (dimensionless)
    Cs_sq = 267.1
    # Degeneracy factor for symmetric nuclear matter
    nu = 4
    # Fermi momentum in fm^-1
    kF_fm = 1.42
    # Radius of heavy nuclei in fm
    R_fm = 5.5
    # Effective mass ratio at saturation point
    M_star_ratio = 0.78
    # Conversion constant hbar*c in MeV*fm
    hbar_c = 197.3

    # --- Derived quantities ---
    # Effective mass in MeV
    M_star = M_star_ratio * M
    # Fermi momentum in MeV
    kF = kF_fm * hbar_c
    # Dimensionless product of Fermi momentum and radius
    kF_R = kF_fm * R_fm

    # --- Calculation ---
    # The term derived from the momentum integral
    kF_over_M_star = kF / M_star
    integral_term_val = math.asinh(kF_over_M_star) - kF / math.sqrt(kF**2 + M_star**2)

    # The full expression for the fractional change per unit of eta
    # (Delta M* / M*) / eta
    # Note: the C_s^2/M^2 part makes the prefactor dimensionless if M is in MeV
    # The formula is (Cs^2/M^2) * (nu * kF*R / (4*pi^2)) * integral_term
    
    prefactor = (Cs_sq / M**2) * (nu * kF_R / (4 * math.pi**2))
    
    relative_change_per_eta = prefactor * integral_term_val
    
    # Percentage change per unit of eta
    percentage_change_per_eta = relative_change_per_eta * 100

    # --- Outputting the results ---
    print("This script calculates the percentage change in the effective mass of the Dirac nucleon per unit of the NC parameter eta.")
    print("The formula used is:")
    print("(% Change / eta) = [Cs^2/M^2 * nu*(kF*R)/(4*pi^2)] * [arsinh(kF/M*) - kF/sqrt(kF^2 + M*^2)] * 100")
    print("\n--- Numerical values used in the calculation ---")
    print(f"Cs^2 = {Cs_sq}")
    print(f"M = {M:.2f} MeV")
    print(f"nu = {nu}")
    print(f"kF*R = {kF_R:.4f}")
    print(f"pi = {math.pi:.4f}")
    print(f"kF = {kF:.4f} MeV")
    print(f"M* = {M_star:.4f} MeV")
    print("\n--- Final Equation with Numbers ---")
    print(f"Percentage Change / eta = [{Cs_sq:.1f}/{M:.2f}^2 * {nu}*({kF_R:.4f})/(4*{math.pi:.4f}^2)] * [arsinh({kF:.4f}/{M_star:.4f}) - {kF:.4f}/sqrt({kF:.4f}^2 + {M_star:.4f}^2)] * 100")
    
    print("\n--- Result ---")
    print(f"The calculated percentage change per unit of eta (in MeV^-2) is: +{percentage_change_per_eta:.6f}%")

if __name__ == '__main__':
    calculate_mass_change()