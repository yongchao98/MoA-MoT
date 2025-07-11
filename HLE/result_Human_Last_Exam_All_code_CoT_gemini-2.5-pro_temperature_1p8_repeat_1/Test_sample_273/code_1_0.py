import numpy as np

def demonstrate_instability():
    """
    Demonstrates the instability of a 3D soliton with only exchange and DMI.

    The total energy of a soliton of size lambda is modeled by:
    E(lambda) = A * lambda - B * lambda^2

    For a stationary point (energy extremum), dE/dlambda = A - 2*B*lambda = 0.
    This occurs at lambda_0 = A / (2*B).

    The second derivative d^2E/dlambda^2 = -2*B. Since B (from DMI) must be
    positive for a soliton to form, the second derivative is negative.
    This proves the stationary point is an energy MAXIMUM, not a stable minimum.

    We will choose A=2.0 and B=1.0 to place this maximum at lambda = 1.0 for clarity.
    """

    # Constants for the energy equation E(lambda) = A*lambda - B*lambda^2
    # These represent the initial exchange and DMI energy contributions.
    A = 2.0  # Proportional to the initial exchange energy
    B = 1.0  # Proportional to the magnitude of the initial DMI energy

    # The normalized size for which the energy is at a maximum
    lambda_max = A / (2 * B)

    # The energy function
    def energy(lam):
        return A * lam - B * lam**2

    # Calculate the energy at the maximum and at nearby points
    energy_at_max = energy(lambda_max)
    
    # A smaller size (collapse tendency)
    lambda_small = 0.5 * lambda_max
    energy_small = energy(lambda_small)

    # A larger size (expansion tendency)
    lambda_large = 2.0 * lambda_max
    energy_large = energy(lambda_large)
    
    # A very large size to show the unbounded nature
    lambda_very_large = 10.0 * lambda_max
    energy_very_large = energy(lambda_very_large)

    print("--- Soliton Instability Demonstration ---")
    print(f"The energy of the soliton is modeled by the equation:")
    print(f"E(lambda) = {A:.2f} * lambda - {B:.2f} * lambda^2\n")

    print(f"The energy has an extremum (a maximum) at a size lambda = {lambda_max:.2f}.")
    print(f"Energy at this maximum E({lambda_max:.2f}) = {energy_at_max:.4f}\n")

    print("Testing stability by changing the size:")
    print(f"If the soliton shrinks to lambda = {lambda_small:.2f}, its energy becomes E = {energy_small:.4f}.")
    print(f"If the soliton expands to lambda = {lambda_large:.2f}, its energy becomes E = {energy_large:.4f}.")
    print(f"If the soliton expands further to lambda = {lambda_very_large:.2f}, its energy becomes E = {energy_very_large:.4f}.\n")

    print("Conclusion: The soliton is unstable.")
    print("Any deviation from the size lambda = 1.0 leads to a lower energy state.")
    print("The energy is unbounded below as lambda increases, predicting catastrophic expansion.")


if __name__ == '__main__':
    demonstrate_instability()