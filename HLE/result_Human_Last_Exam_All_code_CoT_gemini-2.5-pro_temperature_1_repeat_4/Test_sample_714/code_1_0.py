import numpy as np

def calculate_sodium_lifetime():
    """
    Calculates the theoretical lifetime of the Sodium 3p state using a hydrogenic model
    and compares it to the experimental value.
    """
    # Step 1: Define constants and given values from the problem.
    # Physical constants
    c = 2.998e8  # Speed of light in m/s
    hbar = 1.054e-34  # Reduced Planck constant in J*s
    q = 1.602e-19  # Elementary charge in C
    epsilon_0 = 8.854e-12  # Permittivity of free space in F/m
    a0 = 5.29e-11  # Bohr radius in m

    # Given values from the problem
    lambda_m = 589e-9  # Wavelength in m
    tau_exp_s = 16.2e-9  # Experimental lifetime in s
    
    # We use an effective nuclear charge Z=1 for the hydrogenic approximation of the valence electron.
    Z = 1
    
    # The upper state is 3p (l=1), so its degeneracy g_u is 2*l+1.
    l_upper = 1
    g_u = 2 * l_upper + 1 # This is g2 in the problem description

    print("--- Input Values ---")
    print(f"Wavelength (lambda): {lambda_m * 1e9:.1f} nm")
    print(f"Experimental lifetime (tau_exp): {tau_exp_s * 1e9:.1f} ns")
    print(f"Bohr radius (a0): {a0:.2e} m")
    print(f"Effective nuclear charge (Z): {Z}")
    print(f"Degeneracy of upper state (g_u): {g_u}")
    print("-" * 20)

    # Step 2: Calculate the angular frequency (omega) of the transition.
    omega = 2 * np.pi * c / lambda_m

    # Step 3: Calculate the magnitude of the radial integral for the 3p -> 3s transition.
    # The analytical solution for the integral |∫ r^3 * R_3,1 * R_3,0 dr| is 9*sqrt(2)*a0/Z.
    radial_integral_mag = 9 * np.sqrt(2) * a0 / Z

    # Step 4: Calculate the line strength (S).
    # S = max(l_upper, l_lower) * (radial_integral)^2. Here l_lower=0, l_upper=1, so max=1.
    line_strength = (radial_integral_mag)**2

    # Step 5: Calculate the Einstein A coefficient (A_21), the spontaneous emission rate.
    # A_ul = (e^2 * omega^3) / (3 * pi * epsilon_0 * hbar * c^3) * S / g_u
    numerator = q**2 * omega**3 * line_strength
    denominator = 3 * np.pi * epsilon_0 * hbar * c**3 * g_u
    A_21 = numerator / denominator

    # Step 6: Calculate the theoretical lifetime (tau_th), which is the reciprocal of A_21.
    tau_th_s = 1 / A_21
    
    # Step 7: Compare the theoretical lifetime to the experimental lifetime.
    ratio = tau_th_s / tau_exp_s
    
    print("--- Calculation Results ---")
    print(f"Angular frequency (omega): {omega:.3e} rad/s")
    print(f"Radial integral magnitude: {radial_integral_mag:.3e} m")
    print(f"Line strength (S): {line_strength:.3e} m^2")
    print(f"Einstein A coefficient (A_21): {A_21:.3e} s^-1")
    print(f"Calculated theoretical lifetime (tau_th): {tau_th_s * 1e9:.3f} ns")
    print("-" * 20)

    print("--- Final Comparison ---")
    final_eq_str = f"Ratio = tau_th / tau_exp = {tau_th_s * 1e9:.3f} ns / {tau_exp_s * 1e9:.1f} ns = {ratio:.3f}"
    print(final_eq_str)
    
    print("\nThe calculated ratio is approximately 0.1, meaning the theoretical lifetime is one-tenth as long as the experiment.")

calculate_sodium_lifetime()