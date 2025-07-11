import math

def calculate_lindhard_function():
    """
    Calculates the Lindhard polarization function at q=0, omega=0.

    The value is equivalent to the negative of the density of states at the
    Fermi level, g(epsilon_F), for a 3D electron gas. This value depends on
    the electron density 'n' of the gas.

    To provide a specific numerical value, we assume the electron density
    of Sodium (Na), a common example for the free electron gas model.
    """

    # Physical constants in SI units
    hbar = 1.054571817e-34  # Reduced Planck constant (J*s)
    m = 9.10938356e-31      # Electron mass (kg)
    pi = math.pi

    # Assumed parameter: electron density for Sodium (Na) in m^-3
    n = 2.65e28

    # Step 1: Calculate the Fermi momentum k_F
    # k_F = (3 * pi^2 * n)^(1/3)
    k_F = (3 * pi**2 * n)**(1/3)

    # Step 2: Calculate the density of states at the Fermi level g(epsilon_F)
    # g(epsilon_F) = (m * k_F) / (pi^2 * hbar^2)
    g_epsilon_F = (m * k_F) / (pi**2 * hbar**2)

    # Step 3: The Lindhard function at (q=0, omega=0) is -g(epsilon_F)
    lindhard_value = -g_epsilon_F

    # Print the explanation and the final equation with all numbers
    print("The Lindhard function at zero frequency and zero momentum transfer is Π(0,0) = -g(ε_F).")
    print("This value depends on the electron density 'n'. We assume the density of Sodium (n = {:.2e} m^-3).".format(n))
    print("\nThe calculation proceeds as follows:")
    print("Π(0,0) = -g(ε_F) = - (m * k_F) / (π² * ħ²)")
    
    # Substituting k_F = (3 * π² * n)^(1/3)
    print("         = - (m * (3 * π² * n)^(1/3)) / (π² * ħ²)\n")

    # Final equation with numerical values
    print("Final Equation:")
    print("Π(0,0) = - ({:.4e} kg * ({:.4e} m^-1)) / ({:.5f} * ({:.4e} J·s)²)".format(m, k_F, pi**2, hbar))
    print("       = {:.4e} J⁻¹·m⁻³".format(lindhard_value))

# Execute the calculation
calculate_lindhard_function()
<<< -7.6565e+46>>>