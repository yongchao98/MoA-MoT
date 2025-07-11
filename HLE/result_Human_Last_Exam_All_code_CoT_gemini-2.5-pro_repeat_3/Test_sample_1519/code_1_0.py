import math

def calculate_flow_speed():
    """
    Calculates the flow speed at which the pressure at the bottom of a river becomes zero.
    """
    # Given constants
    rho = 1000  # Density of water in kg/m^3
    g = 10      # Acceleration due to gravity in m/s^2
    H = 10      # Depth of the river in meters

    # The problem requires the pressure at the bottom to be zero.
    # This happens when the initial hydrostatic pressure (rho * g * H)
    # is equal to the dynamic pressure drop (1/2 * rho * v^2).
    #
    # Equation:
    # rho * g * H = 1/2 * rho * v^2
    #
    # Simplifying by removing rho from both sides:
    # g * H = 1/2 * v^2
    #
    # Solving for v:
    # v^2 = 2 * g * H
    # v = sqrt(2 * g * H)

    # Calculate the final velocity
    v_squared = 2 * g * H
    v = math.sqrt(v_squared)

    # Print the calculation steps
    print("To find the flow speed 'v' where the pressure at the bottom is zero, we set the hydrostatic pressure equal to the dynamic pressure.")
    print("The equation is: rho * g * H = (1/2) * rho * v^2")
    print("This simplifies to: g * H = (1/2) * v^2")
    print("Solving for v, we get: v = sqrt(2 * g * H)")
    print("\nSubstituting the given values:")
    print(f"v = sqrt(2 * {g} * {H})")
    print(f"v = sqrt({v_squared})")
    print(f"\nThe required flow speed v is: {v} m/s")

calculate_flow_speed()
<<<14.142135623730951>>>