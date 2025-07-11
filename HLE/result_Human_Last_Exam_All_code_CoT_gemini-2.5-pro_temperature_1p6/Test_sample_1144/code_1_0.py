import math

# This script calculates the ratio of the uncertainty of momentum to the momentum
# of an electron in the first Bohr orbit.

def solve_momentum_ratio():
    """
    Calculates and prints the ratio of momentum uncertainty to momentum
    for an electron in the first Bohr orbit.
    """
    # Constants and given values
    # Bohr radius (a0) in meters. It is the radius of the first Bohr orbit.
    bohr_radius = 5.29177210903e-11  # meters

    # Given uncertainty in position (delta_x) in meters.
    # 10 pm = 10 * 10^-12 m
    uncertainty_in_position = 10e-12  # meters

    # According to the Heisenberg Uncertainty Principle, the minimum uncertainty
    # in momentum (delta_p) is given by: delta_p = ħ / (2 * delta_x)
    # where ħ is the reduced Planck constant.

    # According to the Bohr model, the momentum (p) of an electron in the
    # first orbit is given by: p = ħ / a0, where a0 is the Bohr radius.

    # We need to find the ratio: delta_p / p
    # Ratio = (ħ / (2 * delta_x)) / (ħ / a0)
    # The ħ terms cancel out, simplifying the expression to:
    # Ratio = a0 / (2 * delta_x)

    # Calculate the ratio
    ratio = bohr_radius / (2 * uncertainty_in_position)

    # Print the final equation with each number and the result
    print("The problem is to find the ratio: (Uncertainty in Momentum) / (Momentum)")
    print("The simplified final equation is: Bohr Radius / (2 * Uncertainty in Position)")
    print(f"{bohr_radius} / (2 * {uncertainty_in_position}) = {ratio}")

# Execute the function
solve_momentum_ratio()