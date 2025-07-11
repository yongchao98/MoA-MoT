import math

def solve_flow_speed():
    """
    Calculates the flow speed 'v' at which the pressure at the bottom of a river becomes zero.
    """
    # Given values
    rho = 1000  # Density of water in kg/m^3
    g = 10      # Acceleration due to gravity in m/s^2
    H = 10      # Depth of the river in meters

    # The initial gauge pressure at the bottom is P_static = rho * g * H.
    # The pressure drop due to flow is the dynamic pressure, P_dynamic = 1/2 * rho * v^2.
    # We set P_static = P_dynamic to find the speed 'v' where the net pressure is zero.
    # rho * g * H = 1/2 * rho * v^2
    # This simplifies to: g * H = 1/2 * v^2
    # Solving for v: v = sqrt(2 * g * H)

    # Calculate the result
    v = math.sqrt(2 * g * H)

    # Print the explanation and the final equation with numbers
    print("To find the flow speed 'v' where the pressure at the bottom becomes zero, we set the static pressure equal to the dynamic pressure.")
    print("The equation is: rho * g * H = (1/2) * rho * v^2")
    print("This simplifies to: g * H = (1/2) * v^2")
    print("Solving for v, we get: v = sqrt(2 * g * H)")
    print("\nSubstituting the given values:")
    # The final equation with each number outputted
    print(f"v = sqrt(2 * {g} * {H})")
    print(f"v = sqrt({2 * g * H})")
    print(f"\nThe required flow speed v is: {v:.4f} m/s")

solve_flow_speed()