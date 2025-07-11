import math

def calculate_lifetime_and_ratio():
    """
    Calculates the theoretical lifetime of the Sodium 3p state using a hydrogenic model
    and compares it to the experimental value.
    """
    # Constants
    c = 299792458  # Speed of light in m/s
    a0 = 5.29177210903e-11  # Bohr radius in m
    alpha = 7.2973525693e-3  # Fine-structure constant
    
    # Parameters from the problem
    Z = 11  # Nuclear charge of Sodium
    lambda_nm = 589  # Wavelength in nm
    lambda_m = lambda_nm * 1e-9  # Wavelength in m
    tau_exp_ns = 16.2  # Experimental lifetime in ns
    tau_exp_s = tau_exp_ns * 1e-9  # Experimental lifetime in s
    
    # Per problem instruction, assume g2/g1 = 2. Since g1 (for 3s, l=0) is 2*0+1=1, g2=2.
    g2 = 2
    
    # --- Step 1: Calculate angular frequency omega ---
    omega = 2 * math.pi * c / lambda_m
    
    # --- Step 2: Calculate the squared radial integral |I|^2 ---
    # From integrating the provided wavefunctions, the radial integral I = -3*sqrt(2)*a0/Z.
    # The line strength in atomic units of length squared is l_max * |I|^2.
    # For 3p -> 3s, l_max = 1.
    # |I|^2 = (3*sqrt(2)*a0/Z)^2 = 18 * a0^2 / Z^2
    I_sq = 18 * (a0**2) / (Z**2)
    
    # --- Step 3: Calculate the spontaneous emission rate A_21 ---
    # A_21 = (24 * alpha * omega^3 * a0^2) / (c^2 * g2 * Z^2) is not quite right.
    # The correct formula is A_21 = (4 * alpha * omega^3) / (3 * c^2 * g2) * (l_max * |I|^2)
    # A_21 = (4 * alpha * omega^3) / (3 * c^2 * g2) * (1 * 18 * a0^2 / Z^2)
    # A_21 = (24 * alpha * omega^3 * a0^2) / (c^2 * g2 * Z^2)
    
    numerator = 24 * alpha * (omega**3) * (a0**2)
    denominator = (c**2) * g2 * (Z**2)
    A_21 = numerator / denominator
    
    # --- Step 4: Calculate the theoretical lifetime tau_theo ---
    tau_theo_s = 1 / A_21
    tau_theo_ns = tau_theo_s * 1e9
    
    # --- Step 5: Calculate the ratio ---
    ratio = tau_theo_s / tau_exp_s
    
    # --- Print the results ---
    print("--- Calculation Steps ---")
    print(f"1. Angular frequency (ω): {omega:.4e} rad/s")
    print(f"2. Squared radial integral |I|²: {I_sq:.4e} m^2")
    print(f"3. Spontaneous emission rate (A_21): {A_21:.4e} s^-1")
    print("\n--- Final Equation for Theoretical Lifetime (τ_theo) ---")
    print("τ_theo = 1 / A_21")
    print(f"τ_theo = 1 / ( (24 * {alpha:.4e} * ({omega:.4e})³ * ({a0:.4e})²) / (({c:.4e})² * {g2} * {Z}²) )")
    print(f"τ_theo = 1 / ( ({numerator:.4e}) / ({denominator:.4e}) )")
    print(f"τ_theo = 1 / {A_21:.4e} s")
    
    print("\n--- Results ---")
    print(f"Theoretical lifetime (τ_theo): {tau_theo_ns:.1f} ns")
    print(f"Experimental lifetime (τ_exp): {tau_exp_ns:.1f} ns")
    print(f"The theoretical lifetime is {ratio:.1f} times the experimental lifetime.")
    
    # Determine the closest answer choice
    choices = {'A': 30, 'B': 2, 'C': 0.1, 'D': 1, 'E': 100, 'F': 0.5, 'G': 10, 'H': 5}
    closest_choice = min(choices.keys(), key=lambda k: abs(choices[k] - ratio))
    
    print(f"\nThis result is closest to the answer choice of {choices[closest_choice]}.")
    print(f"<<<{closest_choice}>>>")

calculate_lifetime_and_ratio()