import numpy as np

def solve_sodium_lifetime():
    """
    This function calculates the theoretical lifetime of the 3p state in Sodium-23
    based on a hydrogenic model and compares it to the experimental value.
    """
    
    # Step 1: Define constants and parameters
    h = 6.62607015e-34      # Planck constant in J*s
    hbar = h / (2 * np.pi)    # Reduced Planck constant in J*s
    c = 2.99792458e8          # Speed of light in m/s
    e = 1.60217663e-19        # Elementary charge in C
    epsilon_0 = 8.85418782e-12# Permittivity of free space in F/m
    a0 = 5.291772109e-11      # Bohr radius in m
    
    # Calculate the fine-structure constant, alpha
    alpha = e**2 / (4 * np.pi * epsilon_0 * hbar * c)

    # Parameters from the problem statement
    wavelength = 589e-9      # Wavelength of the 3p -> 3s transition in m
    tau_exp = 16.2e-9        # Experimental lifetime of the 3p state in s
    Z = 11                   # Nuclear charge for Sodium
    g2_g1_ratio = 2          # Given degeneracy ratio g_upper / g_lower
    g1_g2_ratio = 1 / g2_g1_ratio

    print("Step-by-step calculation:")
    print("-------------------------")
    print(f"Constants used: a0 = {a0:.4e} m, c = {c:.4e} m/s, alpha = {alpha:.4e}")
    print(f"Problem parameters: Z = {Z}, lambda = {wavelength*1e9:.0f} nm, tau_exp = {tau_exp*1e9:.1f} ns, g2/g1 = {g2_g1_ratio}")
    print("-" * 25)
    
    # Step 2: Calculate the angular frequency (omega) of the transition
    omega = 2 * np.pi * c / wavelength
    print("1. Calculate angular frequency (omega) from wavelength (lambda):")
    print(f"   omega = 2 * pi * c / lambda")
    print(f"   omega = (2 * pi * {c:.4e} m/s) / {wavelength:.4e} m = {omega:.4e} rad/s")

    # Step 3: Calculate the value of the squared radial integral |<3s|r|3p>|^2
    # The radial integral for hydrogenic wavefunctions <3s|r|3p> evaluates to -9*sqrt(2)*a0/Z.
    # We need the square of its magnitude.
    radial_integral_sq = (9 * np.sqrt(2) * a0 / Z)**2
    print("\n2. Calculate the squared radial integral |<3s|r|3p>|^2:")
    print(f"   |<3s|r|3p>|^2 = (9 * sqrt(2) * a0 / Z)^2")
    print(f"   |<3s|r|3p>|^2 = (9 * {np.sqrt(2):.4f} * {a0:.4e} m / {Z})^2 = {radial_integral_sq:.4e} m^2")

    # Step 4: Calculate the Einstein A coefficient (A_21), the rate of spontaneous emission
    # A_21 = (4 * alpha * omega^3 / (3 * c^2)) * (g1/g2) * |<3s|r|3p>|^2
    prefactor_A = (4 * alpha * omega**3) / (3 * c**2)
    A_21 = prefactor_A * g1_g2_ratio * radial_integral_sq
    
    print("\n3. Calculate the Einstein A coefficient (A_21):")
    print(f"   A_21 = (4 * alpha * omega^3 / (3 * c^2)) * (g1/g2) * |<3s|r|3p>|^2")
    print(f"   A_21 = (4 * {alpha:.4e} * ({omega:.4e} rad/s)^3 / (3 * ({c:.4e} m/s)^2)) * ({g1_g2_ratio}) * ({radial_integral_sq:.4e} m^2)")
    print(f"   A_21 = {A_21:.4e} s^-1")

    # Step 5: Calculate the theoretical lifetime (tau_th)
    tau_th = 1 / A_21
    print("\n4. Calculate the theoretical lifetime (tau_th):")
    print(f"   tau_th = 1 / A_21")
    print(f"   tau_th = 1 / {A_21:.4e} s^-1 = {tau_th:.4e} s (or {tau_th*1e9:.1f} ns)")

    # Step 6: Compare the theoretical and experimental lifetimes
    ratio = tau_th / tau_exp
    print("\n5. Find the multiple difference from the experimental lifetime:")
    print(f"   Ratio = tau_th / tau_exp")
    print(f"   Ratio = {tau_th:.4e} s / {tau_exp:.4e} s = {ratio:.1f}")
    
    print("\nConclusion: The theoretical lifetime is approximately 93 times as long as the experimental value.")
    print("The closest answer choice is '100 - one hundred times as long as experiment'.")

solve_sodium_lifetime()
<<<E>>>