import math

def solve_flow_speed():
    """
    Calculates the flow speed 'v' at which the pressure at the bottom of a river becomes zero.
    """
    # Given constants
    rho = 1000  # Density of water in kg/m^3
    g = 10      # Acceleration due to gravity in m/s^2
    H = 10      # Depth of the river in meters

    # The problem requires finding the speed 'v' where the pressure at the bottom is zero.
    # This occurs when the initial hydrostatic pressure (ρ * g * H)
    # is completely converted to dynamic pressure (1/2 * ρ * v^2).
    # The equation is: ρ * g * H = 1/2 * ρ * v^2

    # Calculate the value of v
    # v = sqrt(2 * g * H)
    v = math.sqrt(2 * g * H)

    # Print the explanation and the final equation with numerical values
    print("According to Bernoulli's principle, for the pressure at the bottom to become zero,")
    print("the initial hydrostatic pressure must be equal to the dynamic pressure.")
    print("\nThe equation is:")
    print("ρ * g * H = (1/2) * ρ * v^2")
    print("\nPlugging in the values:")
    print(f"{rho} * {g} * {H} = (1/2) * {rho} * v^2")
    
    # Calculate intermediate values for the equation printout
    static_pressure = rho * g * H
    dynamic_pressure_coeff = 0.5 * rho
    
    print(f"{static_pressure} = {dynamic_pressure_coeff} * v^2")

    print(f"\nSolving for v, we get v = sqrt(2 * g * H)")
    print(f"The required flow speed is: {v:.2f} m/s")

solve_flow_speed()