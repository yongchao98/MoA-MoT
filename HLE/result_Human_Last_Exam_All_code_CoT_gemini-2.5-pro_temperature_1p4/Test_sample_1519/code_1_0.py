import math

def calculate_flow_speed():
    """
    Calculates the flow speed 'v' at which the pressure at the bottom of a river becomes zero.
    """
    # Define the constants based on the problem description
    rho = 1000  # Density of water in kg/m^3
    g = 10      # Acceleration due to gravity in m/s^2, as used in the example
    H = 10      # Depth of the river in meters

    # The problem asks for the flow speed 'v' at which the pressure at the bottom becomes zero.
    # This occurs when the initial static gauge pressure (rho * g * H) is completely
    # converted into dynamic pressure (1/2 * rho * v^2).
    #
    # We set the static pressure equal to the dynamic pressure to find the required speed:
    # rho * g * H = 1/2 * rho * v^2
    #
    # This equation simplifies to:
    # v = sqrt(2 * g * H)

    # Now, we calculate the value of v using the given constants.
    v_squared = 2 * g * H
    v = math.sqrt(v_squared)

    # Print the explanation, the final equation with numbers, and the result.
    print("To find the speed 'v', we set the static pressure equal to the dynamic pressure:")
    print("rho * g * H = 1/2 * rho * v^2")
    print("\nThis simplifies to the equation for 'v':")
    print("v = sqrt(2 * g * H)")
    print("\nPlugging in the given values:")
    print(f"v = sqrt(2 * {g} * {H})")
    print(f"v = sqrt({v_squared})")
    print(f"\nThe required flow speed 'v' is: {v:.4f} m/s")

if __name__ == "__main__":
    calculate_flow_speed()