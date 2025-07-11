import scipy.constants

def calculate_momentum_ratio():
    """
    Calculates the ratio of the uncertainty of the momentum of an electron to its 
    momentum in the first Bohr orbit.

    The ratio is derived from the Heisenberg Uncertainty Principle and the Bohr model,
    and simplifies to: ratio = a0 / (2 * delta_x)
    where:
    - a0 is the Bohr radius.
    - delta_x is the uncertainty in position.
    """

    # Given uncertainty in position in meters
    # 10 pm = 10 * 10^-12 m = 1e-11 m
    delta_x = 10e-12

    # Get the Bohr radius from scipy constants (in meters)
    a0 = scipy.constants.physical_constants["Bohr radius"][0]

    # Calculate the ratio
    ratio = a0 / (2 * delta_x)

    # Print the explanation and the result
    print("This script calculates the ratio of the momentum uncertainty (Δp) to the momentum (p) of an electron in the first Bohr orbit.")
    print("The formula for the ratio is: Δp/p = a₀ / (2 * Δx)")
    print("\nGiven values:")
    print(f"Bohr radius (a₀): {a0} m")
    print(f"Uncertainty in position (Δx): {delta_x} m")
    
    # Print the final equation with all numbers
    print("\nThe final equation is:")
    print(f"ratio = {a0} / (2 * {delta_x})")
    
    print("\nResult:")
    print(f"The ratio of the uncertainty of the momentum to its momentum is: {ratio}")

if __name__ == "__main__":
    calculate_momentum_ratio()
    # The final answer is also printed in the requested format below
    # The calculation is a0 / (2 * delta_x) = (5.29177210903e-11) / (2 * 10e-12)
    # which equals 2.645886054515
    print("<<<2.645886054515>>>")
