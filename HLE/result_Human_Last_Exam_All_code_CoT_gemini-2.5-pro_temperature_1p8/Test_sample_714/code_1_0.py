import math

def calculate_sodium_lifetime():
    """
    Calculates the theoretical lifetime of the Sodium-23 3p state and compares it
    to the experimental value.
    """
    # --- Step 1: Define Constants ---
    # Provided in the problem
    lambda_wl = 589e-9  # Wavelength in meters (589 nm)
    tau_exp = 16.2e-9    # Experimental lifetime in seconds (16.2 ns)
    Z = 11               # Nuclear charge for Sodium
    a0 = 5.29e-11        # Bohr radius in meters
    # Per the problem's hint, we use a degeneracy for the upper state (g2) of 2.
    # This comes from the fine-structure states 3p(3/2) (g2=4) and 3s(1/2) (g1=2), where g2/g1 = 2.
    # The simplest interpretation for the given non-relativistic wavefunctions is to set g2=2.
    g2 = 2

    # Fundamental physical constants
    q = 1.60e-19         # Elementary charge in Coulombs
    epsilon0 = 8.854e-12 # Vacuum permittivity in F/m
    hbar = 1.054e-34     # Reduced Planck constant in J*s
    c = 3.0e8            # Speed of light in m/s

    print("--- Input Parameters ---")
    print(f"Transition Wavelength (lambda): {lambda_wl*1e9:.2f} nm")
    print(f"Experimental Lifetime (tau_exp): {tau_exp*1e9:.2f} ns")
    print(f"Nuclear Charge (Z): {Z}")
    print(f"Upper State Degeneracy (g2): {g2}")
    print("\n")

    # --- Step 2: Calculate the Radial Integral ---
    # The transition is 3p -> 3s (n=3, l=1 to n=3, l'=0)
    # The radial integral <3,0|r|3,1> for a hydrogenic atom is -9*sqrt(2)*a0/Z
    n = 3
    l = 1
    # radial_integral = - (3/2) * n * math.sqrt(n**2 - l**2) * a0 / Z
    radial_integral = -9 * math.sqrt(2) * a0 / Z
    
    print("--- Dipole Moment Calculation ---")
    print(f"Radial integral <R_3s|r|R_3p>: {radial_integral:.4e} m")

    # The line strength S is e^2 * |<r>|^2 * L_supremum. For l->l-1, L_sup=l.
    L_sup = l # which is 1
    # Total sum of squared dipole moments
    sum_d_squared = q**2 * radial_integral**2 * L_sup
    print(f"Sum of squared dipole moments (S): {sum_d_squared:.4e} C^2*m^2")

    # Average squared dipole moment is the sum divided by initial state degeneracy
    avg_d_squared = sum_d_squared / g2
    print(f"Average squared dipole moment (|d_if|^2): {avg_d_squared:.4e} C^2*m^2")
    print("\n")

    # --- Step 3: Calculate the Einstein A coefficient (Transition Rate) ---
    omega = 2 * math.pi * c / lambda_wl
    
    # A_21 = (omega^3 / (3 * pi * epsilon_0 * hbar * c^3)) * avg_d_squared
    prefactor = (omega**3) / (3 * math.pi * epsilon0 * hbar * c**3)
    A21 = prefactor * avg_d_squared
    
    print("--- Lifetime Calculation ---")
    print(f"Angular frequency (omega): {omega:.4e} rad/s")
    print(f"Einstein A coefficient (A21): {A21:.4e} s^-1")
    
    # --- Step 4: Calculate the Theoretical Lifetime ---
    tau_th = 1 / A21
    print(f"Theoretical lifetime (tau_th): {tau_th*1e9:.2f} ns")
    print("\n")

    # --- Step 5: Compare Theoretical and Experimental Lifetimes ---
    multiple = tau_th / tau_exp
    print("--- Final Result ---")
    print(f"The theoretical lifetime is calculated to be {tau_th*1e9:.1f} ns.")
    print(f"The experimental lifetime is given as {tau_exp*1e9:.1f} ns.")
    print("The multiple difference is calculated as:")
    print(f"Multiple = tau_th / tau_exp = {tau_th*1e9:.1f} ns / {tau_exp*1e9:.1f} ns = {multiple:.2f}")

    print("\nThis calculated multiple of {:.2f} is closest to the answer choice of 10.".format(multiple))


calculate_sodium_lifetime()