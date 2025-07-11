import math

def calculate_deflection():
    """
    Calculates the membrane's deflection at x = 0.

    The deflection y(0) is derived by solving the differential equation
    parametrically. The final expression for y(0) is (3/4) * (9/4)^(3/5).
    """

    # Define the components of the final equation
    coefficient = 3/4
    base = 9/4
    exponent = 3/5

    # Calculate the final value
    deflection_at_zero = coefficient * (base ** exponent)

    # Print the equation with its numerical components
    print(f"The calculation for the deflection y(0) is based on the equation:")
    print(f"y(0) = {coefficient} * ({base}) ** ({exponent})")

    # Print the final result
    print("\nThe membrane's deflection at x = 0 is:")
    print(deflection_at_zero)

if __name__ == "__main__":
    calculate_deflection()