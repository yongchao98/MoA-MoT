import cmath

def calculate_bare_greens_function():
    """
    Calculates the bare Green's function G_0 for a given single-particle
    energy eigenvalue epsilon_k and various probe frequencies omega.

    The formula used is G_0(k, omega) = 1 / (omega - epsilon_k + i*eta),
    which corresponds to the propagator for an unoccupied particle state.
    """

    # --- Parameters ---
    # Define the single-particle energy eigenvalue (e.g., in electron-volts, eV).
    epsilon_k = 2.5

    # Define a small positive value for eta. In theory, it is an
    # infinitesimal. We use a small number for demonstration.
    eta = 0.1

    # Define the probe frequencies (omega) to test (e.g., in eV).
    # We choose values below, at, and above epsilon_k to observe the behavior.
    omegas = [1.5, 2.5, 3.5]

    print("Calculating the Bare Green's function G_0(k, omega)")
    print(f"for a single-particle state with energy eigenvalue epsilon_k = {epsilon_k} eV.")
    print(f"We use the formula: G_0 = 1 / (omega - epsilon_k + i*eta)")
    print(f"with a small imaginary part eta = {eta} eV.\n")

    # Loop through each omega and calculate G_0
    for omega in omegas:
        # The denominator is a complex number: (omega - epsilon_k) + i*eta
        denominator = complex(omega - epsilon_k, eta)

        # Calculate the Green's function
        G0 = 1 / denominator

        # Print the equation with the specific numbers used
        print(f"For omega = {omega} eV:")
        real_part = omega - epsilon_k
        print(f"  Equation: G_0 = 1 / ({real_part:.1f} + {eta}i)")

        # Print the final complex result, formatted for clarity
        print(f"  Result:   G_0 = {G0:.4f}\n")

# Run the calculation
calculate_bare_greens_function()