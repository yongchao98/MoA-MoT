import numpy as np

def calculate_absorption_cross_section():
    """
    Calculates and prints the equations for the absorption cross-section
    for a molecular chain under two conditions: non-interacting and
    nearest-neighbor interacting molecules.
    """
    # --- 1. Define physical constants (SI units) ---
    hbar = 1.05457e-34  # J*s
    c = 2.99792e8       # m/s
    epsilon_0 = 8.85418e-12 # F/m
    pi = np.pi

    # --- 2. Define system and laser parameters ---
    N = 10  # Number of molecules in the chain
    # Transition energy corresponding to 500 nm wavelength
    lambda_eg = 500e-9 # m
    omega_eg = 2 * pi * c / lambda_eg  # rad/s
    # Typical molecular transition dipole moment (1 Debye)
    mu_eg = 3.33564e-30 # C*m
    # Gaussian pulse duration (ultrashort)
    tau = 50e-15 # s
    # Nearest-neighbor coupling energy (100 cm^-1 converted to Joules)
    J_wavenumber = 100 # cm^-1
    J = J_wavenumber * (100 * hbar * 2 * pi * c) # J

    # --- 3. Calculate the constant prefactor K and tau^2 ---
    # K = 2*pi*tau / (hbar * c * epsilon_0 * sqrt(pi))
    K_factor = (2 * pi * tau) / (hbar * c * epsilon_0 * np.sqrt(pi))
    tau_squared = tau**2

    print("Derivation of the absorption cross-section sigma(omega):")
    print("-" * 60)
    print("The general equation form is:")
    print("sigma(omega) = Prefactor * omega * exp(-(omega - omega_fi)^2 * tau^2)\n")

    # --- 4. Case (a): No interaction between molecules ---
    print("Case a) No interaction between molecules:")
    # The total squared transition dipole is the incoherent sum
    D2_a = N * mu_eg**2
    # The prefactor for the equation
    prefactor_a = K_factor * D2_a
    # The transition frequency is just that of a single molecule
    omega_fi_a = omega_eg

    print("The absorption is the sum of N identical, degenerate transitions.")
    print("Total squared transition dipole moment = N * |mu_eg|^2")
    print(f"Transition frequency omega_fi = omega_eg = {omega_fi_a:.3e} rad/s")
    print("\nThe final equation is:")
    print(f"sigma_a(omega) = {prefactor_a:.3e} * omega * exp(-(omega - {omega_fi_a:.3e})^2 * {tau_squared:.3e})")
    print("-" * 60)

    # --- 5. Case (b): Nearest-neighbor interaction ---
    print("Case b) Nearest-neighbor interaction:")
    # Interaction creates a band of excitons. Only the k=0 state is optically active.
    # The squared transition dipole for the k=0 state is coherently enhanced
    D2_b = N * mu_eg**2 # Same value, but from |sqrt(N)*mu_eg|^2
    # The prefactor for the equation
    prefactor_b = K_factor * D2_b
    # The frequency of the k=0 exciton is shifted by the coupling J
    omega_fi_b = omega_eg + (2 * J / hbar)

    print("Interaction creates a band of exciton states.")
    print("A selection rule allows transitions only to the k=0 exciton state.")
    print("Squared transition dipole moment for this state = |sqrt(N)*mu_eg|^2 = N * |mu_eg|^2 (Superradiance)")
    print(f"Transition frequency omega_fi = omega_eg + 2J/hbar = {omega_fi_b:.3e} rad/s")
    print("\nThe final equation is:")
    print(f"sigma_b(omega) = {prefactor_b:.3e} * omega * exp(-(omega - {omega_fi_b:.3e})^2 * {tau_squared:.3e})")
    print("-" * 60)

calculate_absorption_cross_section()