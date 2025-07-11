import numpy as np

def calculate_sodium_lifetime():
    """
    Calculates the theoretical lifetime of the Sodium-23 3p state
    using a hydrogenic model and compares it to the experimental value.
    """
    # --- Physical Constants ---
    a0 = 5.292e-11      # Bohr radius in meters (m)
    c = 2.998e8         # Speed of light in m/s
    alpha = 1 / 137.036   # Fine-structure constant

    # --- Given Experimental Parameters ---
    lambda_nm = 589
    lambda_m = lambda_nm * 1e-9  # Wavelength in meters
    tau_exp_ns = 16.2
    tau_exp_s = tau_exp_ns * 1e-9  # Experimental lifetime in seconds

    # --- Model Parameters ---
    # For a hydrogenic model of an alkali atom, the valence electron is outside
    # a screened core. The simplest model assumes perfect screening,
    # so we use an effective nuclear charge Z = 1.
    Z = 1.0
    n = 3
    l_i = 1  # Initial state is 3p (l=1)

    print("Step 1: Calculate transition properties")
    print("---------------------------------------")
    # Angular frequency omega = 2 * pi * c / lambda
    omega = 2 * np.pi * c / lambda_m
    print(f"Angular frequency ω = 2 * π * {c:.3e} / {lambda_m:.3e} = {omega:.3e} rad/s")

    print("\nStep 2: Calculate the radial matrix element <r>")
    print("---------------------------------------------")
    # Using the formula: <r> = (3/2) * n * sqrt(n^2 - l_i^2) * a0 / Z
    I_rad = (3/2) * n * np.sqrt(n**2 - l_i**2) * a0 / Z
    print(f"<r> = (3/2) * {n} * sqrt({n}² - {l_i}²) * {a0:.3e} / {Z} = {I_rad:.3e} m")

    print("\nStep 3: Calculate the theoretical decay rate (Einstein A coefficient)")
    print("-------------------------------------------------------------------")
    # The decay rate is A = (4 * alpha * omega^3) / (9 * c^2) * <r>^2
    # This formula calculates the lifetime for a p -> s transition in a hydrogenic model.
    A_coeff = (4 * alpha * omega**3 * I_rad**2) / (9 * c**2)
    print(f"A = (4 * {alpha:.6f} * ({omega:.3e})³ * ({I_rad:.3e})²) / (9 * ({c:.3e})²) = {A_coeff:.3e} s⁻¹")
    
    print("\nStep 4: Calculate the theoretical lifetime τ_th")
    print("---------------------------------------------")
    # Theoretical lifetime tau_th = 1 / A
    tau_th_s = 1 / A_coeff
    tau_th_ns = tau_th_s * 1e9
    print(f"τ_th = 1 / {A_coeff:.3e} = {tau_th_s:.3e} s = {tau_th_ns:.3f} ns")
    
    print("\nStep 5: Compare theoretical and experimental lifetimes")
    print("-----------------------------------------------------")
    # Ratio = tau_th / tau_exp
    ratio = tau_th_ns / tau_exp_ns
    print(f"Experimental lifetime τ_exp = {tau_exp_ns} ns")
    print(f"The ratio τ_th / τ_exp = {tau_th_ns:.3f} ns / {tau_exp_ns} ns = {ratio:.3f}")
    
    print(f"\nThe theoretical lifetime ({tau_th_ns:.2f} ns) is approximately {ratio:.1f} times as long as the experimental lifetime ({tau_exp_ns} ns).")
    print("This corresponds to 'one tenth as long as experiment'.")

calculate_sodium_lifetime()
<<<C>>>