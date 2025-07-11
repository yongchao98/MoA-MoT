import math

def calculate_momentum_ratio():
    """
    Calculates the ratio of the uncertainty of the momentum of an electron
    to its momentum in the first Bohr orbit.
    """
    # Constants and given values
    # Bohr radius (a₀) in meters
    a0 = 5.29177210903e-11
    # Given uncertainty in position in picometers
    delta_x_pm = 10
    # Convert pm to meters (1 pm = 1e-12 m)
    delta_x_m = delta_x_pm * 1e-12

    # The ratio Δp/p simplifies to a₀ / (2 * Δx)
    ratio = a0 / (2 * delta_x_m)

    print("This script calculates the ratio of momentum uncertainty to momentum (Δp/p) for an electron in the first Bohr orbit.")
    print("The simplified formula used is: a₀ / (2 * Δx)\n")
    print(f"Values used:")
    print(f"Bohr radius (a₀) = {a0:.5e} m")
    print(f"Uncertainty in position (Δx) = {delta_x_pm} pm = {delta_x_m:.5e} m\n")

    print("Final Equation with values:")
    # Using f-string formatting to display the equation with numbers
    print(f"Ratio = {a0} / (2 * {delta_x_m})")

    print(f"\nResult:")
    print(f"The ratio Δp / p is: {ratio}")

# Execute the function
calculate_momentum_ratio()