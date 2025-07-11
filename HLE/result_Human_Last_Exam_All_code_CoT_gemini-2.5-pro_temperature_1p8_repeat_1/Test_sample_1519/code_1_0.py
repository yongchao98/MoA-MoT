import math

def solve_flow_speed():
    """
    Calculates the flow speed 'v' at which the gauge pressure at the
    bottom of a river becomes zero, based on Bernoulli's principle.
    """
    # Given constants
    rho = 1000  # Density of water in kg/m^3
    g = 10      # Acceleration due to gravity in m/s^2
    H = 10      # Depth of the river in meters

    # According to Bernoulli's principle, the initial hydrostatic pressure
    # is converted to dynamic pressure as the water flows.
    # We want to find the velocity 'v' where the gauge pressure at the bottom is zero.
    # This occurs when the dynamic pressure equals the initial hydrostatic pressure.
    #
    # The equation is: Hydrostatic Pressure = Dynamic Pressure
    # ρ * g * H = 1/2 * ρ * v^2
    
    # We solve for v. The density 'ρ' cancels out from both sides.
    # v = sqrt(2 * g * H)
    v = math.sqrt(2 * g * H)

    # Print the final equation with all the numbers, as requested.
    print("The final equation to solve is when hydrostatic pressure equals dynamic pressure:")
    print(f"{rho} * {g} * {H} = 0.5 * {rho} * v^2")

    # Print the solution
    print(f"\nSolving for v, we get v = sqrt(2 * g * H).")
    print(f"The required flow speed is {v:.2f} m/s.")

solve_flow_speed()