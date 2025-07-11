import numpy as np

def calculate_mass_change():
    """
    Calculates the percentage change in the effective mass of the Dirac nucleon
    due to noncommutative geometry.
    """
    # Step 1: Define constants and convert units
    M = 938.93         # Nucleon mass in MeV
    Cs2 = 267.1        # Scalar coupling constant (dimensionless in this model)
    kF_fm_inv = 1.42   # Fermi momentum in fm^-1
    R_fm = 5.5         # Radius of the nucleus in fm
    nu = 4             # Degeneracy factor for nuclear matter
    M_star_factor = 0.78 # M* = 0.78 * M
    hc = 197.3         # Conversion factor h-bar * c in MeV fm

    # Convert units to MeV
    kF = kF_fm_inv * hc
    R = R_fm / hc

    print("--- Step 1: Constants and Parameters (in Natural Units) ---")
    print(f"Nucleon Mass M = {M:.2f} MeV")
    print(f"Scalar Coupling C_s^2 = {Cs2}")
    print(f"Fermi Momentum k_F = {kF_fm_inv} fm^-1 = {kF:.2f} MeV")
    print(f"Nuclear Radius R = {R_fm} fm = {R:.6f} MeV^-1")
    print(f"Degeneracy Factor nu = {nu}")
    print("\n" + "="*50 + "\n")

    # Step 2: Calculate intermediate values
    M_star = M_star_factor * M
    
    kF2 = kF**2
    M_star2 = M_star**2
    
    # Evaluate the definite momentum integral from 0 to k_F
    sqrt_term = np.sqrt(kF2 + M_star2)
    val_at_kF = (kF2 + 2 * M_star2) / sqrt_term
    val_at_0 = (0 + 2 * M_star2) / np.sqrt(0 + M_star2) # This is 2*M_star
    
    # This is the result of the 1D integral part
    P = val_at_kF - val_at_0

    print("--- Step 2: Intermediate Calculations ---")
    print(f"Effective Mass M* = {M_star_factor} * M = {M_star:.2f} MeV")
    print(f"Momentum Integral Result P = {P:.4f} MeV")
    print("\n" + "="*50 + "\n")

    # Step 3: Calculate the mass shift per unit eta
    # Formula: d(M*)/d(eta) = (Cs2 / M^2) * nu * M* * P * (3 * R / (32 * pi^2))
    mass_shift_per_eta = (Cs2 / M**2) * nu * M_star * (3 * R / (32 * np.pi**2)) * P

    print("--- Step 3: Mass Shift per unit of NC parameter eta ---")
    print(f"d(M*)/d(eta) = {mass_shift_per_eta:.6f} MeV^-1")
    print("\n" + "="*50 + "\n")

    # Step 4: Determine a physically motivated value for eta
    # We assume eta is the average momentum squared of a nucleon in the Fermi gas.
    # eta = <p^2> = (3/5) * kF^2
    eta = (3.0/5.0) * kF**2
    
    print("--- Step 4: Estimating the NC parameter eta ---")
    print(f"Assuming eta is the average momentum squared: <p^2> = (3/5) * kF^2")
    print(f"eta = (3/5) * ({kF:.2f} MeV)^2 = {eta:.2f} MeV^2")
    print("\n" + "="*50 + "\n")

    # Step 5: Calculate the final percentage change
    # Delta_M* = eta * (d(M*)/d(eta))
    delta_M_star = eta * mass_shift_per_eta
    
    # Percentage Change = 100 * Delta_M* / M*
    percent_change = 100 * delta_M_star / M_star
    
    print("--- Step 5: Final Percentage Change Calculation ---")
    print("Final Equation:")
    print(f"Percentage Change = 100 * (Delta M*) / M*")
    print(f"Delta M* = eta * (d(M*)/d(eta)) = {eta:.2f} MeV^2 * {mass_shift_per_eta:.6f} MeV^-1 = {delta_M_star:.4f} MeV")
    print(f"M* = {M_star:.2f} MeV")
    print("\nResult:")
    print(f"Percentage Change = 100 * {delta_M_star:.4f} / {M_star:.2f} = {percent_change:.2f}%")

calculate_mass_change()