import cmath

def calculate_bare_greens_function():
    """
    Calculates and prints the bare Green's function G_0 for a given
    frequency, single-particle energy, and a small imaginary part.

    This example assumes a "particle" state, where the energy level is
    unoccupied, so the sign of the imaginary part is positive.
    """
    # Define the parameters for the calculation
    # omega (ω): The frequency (energy) of the propagating particle.
    # epsilon_k (ε_k): The single-particle energy eigenvalue.
    # delta (δ): A small positive infinitesimal to enforce causality.
    omega = 2.5  # units of energy
    epsilon_k = 2.2  # units of energy
    delta = 0.1  # a small positive value

    # The denominator of the Green's function
    denominator = (omega - epsilon_k) + 1j * delta

    # The bare Green's function G_0 is the reciprocal of the denominator
    G0 = 1 / denominator

    # Print the full equation and the result, showing each number used.
    # Greek letters ω, ε, δ are used for clarity.
    print(f"The formula for the bare Green's function is:")
    print(f"G_0(k, \u03C9) = 1 / (\u03C9 - \u03B5_k + i\u03B4)\n")

    print(f"Substituting the values:")
    print(f"G_0 = 1 / ({omega} - {epsilon_k} + i*{delta})")
    print(f"G_0 = 1 / ({denominator.real:.2f} + {denominator.imag:.2f}i)")
    print(f"\nThe final result is:")
    print(f"G_0 = {G0.real:.4f} + ({G0.imag:.4f})i")


# Run the calculation and print the results
calculate_bare_greens_function()