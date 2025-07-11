import math

def solve_flow_speed():
    """
    Calculates the flow speed 'v' at which the pressure at the bottom of a river becomes zero.
    """
    # Given parameters
    rho = 1000  # Density of water in kg/m^3
    g = 10      # Acceleration due to gravity in m/s^2
    H = 10      # Depth of the river in meters

    # The problem states that the initial pressure at the bottom (rho * g * H)
    # is reduced to zero due to the flow. This means the initial pressure
    # is converted entirely into dynamic pressure (1/2 * rho * v^2).
    #
    # Equation: rho * g * H = 1/2 * rho * v^2
    # Simplifying by dividing by rho: g * H = 1/2 * v^2
    # Solving for v: v = sqrt(2 * g * H)

    # Calculate the speed v
    v = math.sqrt(2 * g * H)

    # Print the explanation and the final equation with substituted values
    print("The relationship derived from Bernoulli's principle is: v = sqrt(2 * g * H)")
    print("Substituting the given values:")
    print(f"v = sqrt(2 * {g} * {H})")
    print(f"The calculated flow speed is: {v} m/s")

solve_flow_speed()
print("<<<14.142135623730951>>>")