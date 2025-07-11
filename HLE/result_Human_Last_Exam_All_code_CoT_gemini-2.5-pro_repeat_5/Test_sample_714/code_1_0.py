import numpy as np

def calculate_sodium_lifetime():
    """
    Calculates the theoretical lifetime of the Sodium-23 3p state using a hydrogenic model
    and compares it to the experimental value.
    """
    # --- Constants ---
    q = 1.602176634e-19      # Elementary charge in C
    epsilon_0 = 8.8541878128e-12 # Permittivity of free space in F/m
    hbar = 1.054571817e-34       # Reduced Planck constant in J s
    c = 2.99792458e8            # Speed of light in m/s
    a_0 = 5.29177210903e-11     # Bohr radius in m

    # --- Problem Parameters ---
    lambda_nm = 589              # Wavelength in nm
    lambda_m = lambda_nm * 1e-9  # Wavelength in m
    Z = 11                       # Nuclear charge for the hydrogenic model
    g2_over_g1 = 2.0             # Given degeneracy ratio
    g1_over_g2 = 1.0 / g2_over_g1
    tau_exp_ns = 16.2            # Experimental lifetime in ns
    tau_exp_s = tau_exp_ns * 1e-9 # Experimental lifetime in s

    # --- Transition Properties (3p -> 3s) ---
    # Upper state (2): n=3, l=1. Lower state (1): n=3, l=0.
    l_upper = 1
    l_lower = 0
    l_max = max(l_upper, l_lower)

    # Step 1: Calculate the squared radial integral (I_r)^2
    # The analytical result for the integral I_r = integral(R_3,1 * r * R_3,0 * r^2 dr)
    # using the provided wavefunctions is I_r = -27 * sqrt(2) * a_0 / Z.
    I_r_sq = (27 * np.sqrt(2) * a_0 / Z)**2

    # Step 2: Calculate the angular frequency omega
    omega = 2 * np.pi * c / lambda_m

    # Step 3: Calculate the Einstein A coefficient (A_21)
    # A_21 = (g1/g2) * (q^2 * omega^3 * l_max * I_r^2) / (3 * pi * epsilon_0 * hbar * c^3)
    term1 = g1_over_g2
    term2 = (q**2 * omega**3 * l_max * I_r_sq)
    term3 = (3 * np.pi * epsilon_0 * hbar * c**3)
    A_21 = term1 * term2 / term3

    # Step 4: Calculate the theoretical lifetime (tau_th)
    tau_th_s = 1 / A_21

    # Step 5: Compare theoretical and experimental lifetimes
    ratio = tau_th_s / tau_exp_s

    # --- Final Output ---
    print("The final equation for the theoretical lifetime is tau_th = 1 / A_21, where:")
    print("A_21 = (g1/g2) * (q^2 * omega^3 * l_max * I_r^2) / (3 * pi * epsilon_0 * hbar * c^3)\n")
    print("Using the following values:")
    print(f"  g1/g2 = {g1_over_g2}")
    print(f"  q = {q:.5e} C (electron charge)")
    print(f"  omega = {omega:.5e} rad/s (from lambda = {lambda_nm} nm)")
    print(f"  l_max = {l_max}")
    print(f"  I_r^2 = {I_r_sq:.5e} m^2 (from a_0={a_0:.5e} m, Z={Z})")
    print(f"  pi = {np.pi:.5f}")
    print(f"  epsilon_0 = {epsilon_0:.5e} F/m")
    print(f"  hbar = {hbar:.5e} J.s")
    print(f"  c = {c:.5e} m/s\n")

    print(f"The calculated Einstein coefficient is A_21 = {A_21:.4e} s^-1.")
    print(f"The resulting theoretical lifetime is tau_th = {tau_th_s * 1e9:.2f} ns.")
    print(f"The given experimental lifetime is tau_exp = {tau_exp_ns:.2f} ns.\n")
    print(f"The ratio of theoretical to experimental lifetime is: {ratio:.2f}")
    print(f"The calculated lifetime is approximately {round(ratio)} times the experimental lifetime.")

calculate_sodium_lifetime()
<<<E>>>