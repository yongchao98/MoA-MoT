import math

def calculate_flow_speed():
    """
    Calculates the flow speed at which the pressure at the bottom of a river becomes zero.
    """
    # Given values
    rho = 1000  # Density of water in kg/m^3
    g = 10      # Acceleration due to gravity in m/s^2
    H = 10      # Depth of the river in meters

    # The problem requires the pressure at the bottom to decrease to zero.
    # This means the initial hydrostatic pressure (P = ρgh) must be equal to
    # the dynamic pressure created by the flow (P_dynamic = 1/2 * ρ * v^2).
    # So, we set the two equations equal to each other.
    # ρ * g * H = 1/2 * ρ * v^2

    # We can solve for v. Notice that ρ cancels from both sides:
    # g * H = 1/2 * v^2
    # v^2 = 2 * g * H
    v_squared = 2 * g * H
    v = math.sqrt(v_squared)

    # Output the explanation and the equation with numbers as requested.
    print("To find the speed 'v' where the pressure becomes zero, we set the initial hydrostatic pressure equal to the dynamic pressure.")
    print("The equation is: rho * g * H = 1/2 * rho * v^2")
    print("\nSubstituting the given values:")
    # Using float(rho) to ensure floating point division if needed, and clean output for 0.5
    print(f"{float(rho)} * {float(g)} * {float(H)} = 0.5 * {float(rho)} * v^2")

    # Show the calculation steps
    hydrostatic_pressure = rho * g * H
    rho_half = 0.5 * rho
    print(f"{hydrostatic_pressure} = {rho_half} * v^2")
    print(f"v^2 = {hydrostatic_pressure} / {rho_half}")
    print(f"v^2 = {v_squared}")
    print(f"v = sqrt({v_squared})")
    
    # Print the final answer
    print(f"\nThe required flow speed is: {v} m/s")


calculate_flow_speed()
# The final answer is the numerical value of v
<<<14.142135623730951>>>