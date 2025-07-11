import math

def solve_flow_speed():
    """
    Calculates the flow speed 'v' at which the pressure at the bottom of a river becomes zero.
    """
    # Given values
    rho = 1000  # Density of water in kg/m^3
    g = 10      # Acceleration due to gravity in m/s^2
    H = 10      # Depth of the river in meters

    # Initial static pressure at the bottom
    P_static = rho * g * H

    # The condition is that the pressure at the bottom becomes zero.
    # This means the initial static pressure is completely converted to dynamic pressure.
    # P_static = P_dynamic
    # rho * g * H = (1/2) * rho * v^2

    # Solving for v:
    # v = sqrt(2 * g * H)
    v = math.sqrt(2 * g * H)

    print("The initial static pressure at the bottom of the river is P_static = rho * g * H.")
    print(f"P_static = {rho} * {g} * {H} = {P_static} N/m^2")
    print("\nWhen the water flows, the dynamic pressure is P_dynamic = 1/2 * rho * v^2.")
    print("For the pressure at the bottom to decrease to zero, the dynamic pressure must equal the initial static pressure.")
    
    print("\nWe solve the following equation for v:")
    print(f"{P_static} = 0.5 * {rho} * v^2")
    
    print("\nRearranging the equation to find v:")
    print(f"v^2 = 2 * {g} * {H}")
    print(f"v^2 = {2 * g * H}")
    print(f"v = sqrt({2 * g * H})")
    
    print("\nThe required flow speed is:")
    print(f"v = {v} m/s")

solve_flow_speed()
