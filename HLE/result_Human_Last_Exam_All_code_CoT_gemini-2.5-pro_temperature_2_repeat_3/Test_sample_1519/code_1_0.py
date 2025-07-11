import math

def calculate_flow_speed():
    """
    Calculates the flow speed 'v' at which the pressure at the bottom of a river becomes zero.
    """
    # Given values
    g = 10  # Acceleration due to gravity in m/s^2
    H = 10  # Depth of the river in meters
    rho = 1000 # Density of water in kg/m^3 (Note: this value cancels out in the final formula)

    # The static pressure at the bottom is P_static = rho * g * H.
    # When flowing, this is balanced by the dynamic pressure and the remaining pressure:
    # rho * g * H = P_dynamic + 1/2 * rho * v^2
    # We want to find v when P_dynamic = 0.
    # The equation simplifies to v = sqrt(2 * g * H).
    
    # Calculate the speed v
    v_squared = 2 * g * H
    v = math.sqrt(v_squared)

    print("The relationship between static pressure and dynamic pressure gives the equation:")
    print("rho * g * H = 1/2 * rho * v^2")
    print("\nThis simplifies to the formula for the flow speed v:")
    print("v = sqrt(2 * g * H)")
    print("\nSubstituting the given values:")
    print(f"v = sqrt(2 * {g} * {H})")
    print(f"v = sqrt({v_squared})")
    print(f"\nThe required flow speed is {v} m/s.")

calculate_flow_speed()