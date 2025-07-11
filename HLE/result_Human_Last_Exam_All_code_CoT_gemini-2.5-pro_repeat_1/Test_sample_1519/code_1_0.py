import math

def solve_flow_speed():
    """
    Calculates the flow speed 'v' at which the pressure at the bottom of a river becomes zero.
    """
    # Given values
    rho = 1000  # density of water in kg/m^3
    g = 10      # acceleration due to gravity in m/s^2
    H = 10      # depth of the river in meters

    # The initial hydrostatic pressure at the bottom is P = rho * g * H.
    # The dynamic pressure due to flow is (1/2) * rho * v^2.
    # We want to find the speed 'v' where the pressure at the bottom is zero.
    # This occurs when the dynamic pressure equals the initial hydrostatic pressure.
    
    print("The condition is that the hydrostatic pressure equals the dynamic pressure.")
    print("Equation: rho * g * H = (1/2) * rho * v^2")
    print("\nSince 'rho' is on both sides, we can simplify the equation to:")
    print("g * H = 0.5 * v^2")
    
    print("\nNow, we plug in the values for g and H:")
    # We print the full equation with the numbers for clarity as requested.
    print(f"{g} * {H} = 0.5 * v^2")
    
    # Perform the calculation step-by-step
    gh = g * H
    print(f"{gh} = 0.5 * v^2")
    
    v_squared = gh / 0.5
    print(f"v^2 = {gh} / 0.5")
    print(f"v^2 = {v_squared}")
    
    v = math.sqrt(v_squared)
    print(f"v = sqrt({v_squared})")
    
    print(f"\nThe required flow speed v is {v:.2f} m/s.")

solve_flow_speed()

# The final calculated value
v = math.sqrt(2 * 10 * 10)
print(f'<<<{v:.2f}>>>')