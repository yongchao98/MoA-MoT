import numpy as np

def calculate_decay_angle():
    """
    Calculates the angle between a decay product C and its parent particle A.

    The calculation uses Lorentz transformations on the four-momentum of C from
    the rest frame of A to the lab frame.
    """
    # Given velocity of particle A
    beta_A = 0.95

    # Calculate the Lorentz factor, gamma
    gamma_A = 1 / np.sqrt(1 - beta_A**2)

    # Calculate the term in the denominator of the tan(theta) expression
    denominator_term = np.sqrt(2) * beta_A + 1

    # Calculate the tangent of the angle theta_C
    tan_theta_C = 1 / (gamma_A * denominator_term)

    # Calculate the angle in degrees
    theta_C_deg = np.rad2deg(np.arctan(tan_theta_C))

    # Print the equation and the values used
    print("The angle theta_C is found using the formula:")
    print("theta_C = arctan(1 / (gamma_A * (sqrt(2) * beta_A + 1)))")
    print("\nStep-by-step calculation:")
    print(f"1. The velocity of particle A is beta_A = {beta_A}")
    print(f"2. The Lorentz factor is gamma_A = 1 / sqrt(1 - {beta_A}^2) = {gamma_A:.5f}")
    print(f"3. The denominator term is (sqrt(2) * {beta_A} + 1) = {denominator_term:.5f}")
    print(f"4. The tangent of the angle is tan(theta_C) = 1 / ({gamma_A:.5f} * {denominator_term:.5f}) = {tan_theta_C:.5f}")
    print(f"5. The angle in degrees is theta_C = arctan({tan_theta_C:.5f})")
    print(f"\nFinal Answer: The angle is {theta_C_deg:.3f} degrees.")


if __name__ == "__main__":
    calculate_decay_angle()
