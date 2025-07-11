import math

def calculate_nc_mass_change():
    """
    Calculates the percentage change in the effective mass of the Dirac nucleon
    due to non-commutative geometry.
    """
    # Step 1: Define constants and convert units
    M = 938.93  # Nucleon mass in MeV
    Cs_sq = 267.1  # Scalar coupling constant (dimensionless)
    kF_inv_fm = 1.42  # Fermi momentum in fm^-1
    R_fm = 5.5  # Nuclear radius in fm
    M_star_ratio = 0.78  # Effective mass ratio
    nu = 4  # Degeneracy factor
    h_bar_c = 197.3  # hbar*c in MeV fm

    # Convert to consistent MeV units
    M_star = M_star_ratio * M
    kF = kF_inv_fm * h_bar_c
    R = R_fm / h_bar_c  # R in MeV^-1

    print("--- Input Parameters (in MeV units) ---")
    print(f"Nucleon mass M = {M:.2f} MeV")
    print(f"Effective nucleon mass M* = {M_star:.2f} MeV")
    print(f"Fermi momentum k_F = {kF:.2f} MeV")
    print(f"Nuclear radius R = {R:.6f} MeV^-1")
    print("-" * 20)

    # Step 2: Calculate the term from the integral
    integral_term_num = kF**2 + 2 * M_star**2
    integral_term_den = math.sqrt(kF**2 + M_star**2)
    integral_term_sub = 2 * M_star
    integral_term = integral_term_num / integral_term_den - integral_term_sub
    
    print("--- Intermediate Calculation ---")
    print(f"Integral term [k_F^2 + 2M*^2]/sqrt(k_F^2 + M*^2) - 2M* = {integral_term:.4f} MeV")

    # Step 3: Calculate the fractional change coefficient C
    # C = (M*_NC - M*)/(M* * eta)
    prefactor_C = (Cs_sq / M**2) * (3 * nu * R / (32 * math.pi**2))
    C = prefactor_C * integral_term

    print(f"Fractional change coefficient C = {C:.6g} MeV^-2")
    
    # Step 4 & 5: Define scaling factor S and calculate final percentage change
    # S = (1/R)^2
    S = (1/R)**2
    
    percentage_change = C * S * 100
    
    print("\n--- Final Calculation ---")
    print("The final result is obtained by scaling C with S = (1/R)^2.")
    print(f"Scaling factor S = (1 / {R:.6f})^2 = {S:.2f} MeV^2")
    print(f"Percentage Change = C * S * 100")
    print(f"Percentage Change = {C:.6g} * {S:.2f} * 100 = {percentage_change:.4f}%")
    
    # Compare with answer choices
    print("\nThis result is approximately +0.14%, which is closest to answer choice A.")


calculate_nc_mass_change()
<<<A>>>