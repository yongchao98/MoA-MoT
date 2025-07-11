import numpy as np

def solve_lifetime():
    """
    Calculates the theoretical lifetime of the Sodium 3p state and compares it to the experimental value.
    """
    # Step 1: Define constants and parameters
    q = 1.602e-19  # Electron charge in C
    a0 = 5.29e-11  # Bohr radius in m
    lda = 589e-9   # Wavelength of transition in m
    tau_exp = 16.2e-9 # Experimental lifetime in s
    Z = 1.0        # Effective nuclear charge for the valence electron of Sodium
    c = 3.0e8      # Speed of light in m/s
    h_bar = 1.054e-34 # Reduced Planck constant in J*s
    eps0 = 8.854e-12 # Permittivity of free space in F/m

    # Transition quantum numbers
    n2, l2 = 3, 1  # Upper state: 3p
    n1, l1 = 3, 0  # Lower state: 3s
    l_max = max(l1, l2)

    # Step 2: Calculate angular frequency
    omega = 2 * np.pi * c / lda

    # Step 3: Calculate the radial integral using the analytical result from the provided table
    # I_rad = integral(R_31 * r * R_30 * r^2 dr) from 0 to infinity
    # Based on the provided wavefunctions, this evaluates to: (39*sqrt(2)/4) * a0 / Z
    I_rad = (39 * np.sqrt(2) / 4) * a0 / Z

    # Step 4: Calculate the angular factor for the transition rate
    angular_factor = l_max / (2 * l_max + 1)

    # Step 5: Calculate the Einstein A coefficient (A21)
    A21_prefactor_val = (omega**3 * q**2) / (3 * np.pi * eps0 * h_bar * c**3)
    A21 = A21_prefactor_val * angular_factor * I_rad**2
    
    # Step 6: Calculate the theoretical lifetime
    tau_th = 1 / A21
    
    # Step 7: Calculate the ratio
    ratio = tau_th / tau_exp

    # Output the results, showing the final equation with all numbers
    print("--- Calculation of Theoretical Lifetime (τ_th) ---")
    print("The theoretical lifetime τ_th is the inverse of the Einstein A coefficient (A₂₁).")
    print("A₂₁ = [ (ω³ q²) / (3πε₀ħc³) ] * [ l_max / (2l_max+1) ] * |I_rad|²\n")
    
    print("--- Values Used ---")
    print(f"ω (angular frequency) = {omega:.4e} rad/s")
    print(f"q (electron charge)   = {q:.4e} C")
    print(f"ε₀ (permittivity)     = {eps0:.4e} F/m")
    print(f"ħ (reduced Planck)    = {h_bar:.4e} J·s")
    print(f"c (speed of light)    = {c:.4e} m/s")
    print(f"l_max (angular q.n.)  = {l_max}")
    print(f"I_rad (radial integral) = {I_rad:.4e} m\n")

    print("--- Final Equation for A₂₁ ---")
    eq_str = f"A₂₁ = [ ({omega:.4e}³ * {q:.4e}²) / (3π * {eps0:.4e} * {h_bar:.4e} * {c:.4e}³) ] * [ {l_max} / (2*{l_max}+1) ] * ({I_rad:.4e})²"
    print(eq_str)
    print(f"A₂₁ = {A21:.4e} s⁻¹\n")

    print("--- Final Lifetime Calculation ---")
    print(f"τ_th = 1 / A₂₁ = 1 / {A21:.4e} s⁻¹ = {tau_th:.4e} s")
    print(f"τ_th = {tau_th * 1e9:.2f} ns\n")
    
    print("--- Comparison with Experiment ---")
    print(f"Experimental lifetime τ_exp = {tau_exp * 1e9:.2f} ns")
    print(f"Ratio (τ_th / τ_exp) = {tau_th * 1e9:.2f} ns / {tau_exp * 1e9:.2f} ns = {ratio:.3f}")

    print("\nThe theoretical lifetime is approximately one tenth as long as the experimentally measured lifetime.")

solve_lifetime()
<<<C>>>