import math

def solve_flow_speed():
    """
    Calculates the flow speed at which the pressure at the bottom of a river becomes zero.
    """
    # Step 1: Define the given constants.
    # ρ (rho) is the density of water in kg/m^3.
    # g is the acceleration due to gravity in m/s^2.
    # H is the depth of the river in meters.
    rho = 1000
    g = 10
    H = 10

    # Step 2: Explain the physical principle and the equation to be solved.
    # The pressure at the bottom becomes zero when the initial static pressure (P_static)
    # is completely offset by the dynamic pressure (P_dynamic) from the flow.
    print("To solve this, we set the static pressure equal to the dynamic pressure:")
    print("P_static = P_dynamic")
    print("ρ * g * H = (1/2) * ρ * v^2")
    print("This can be written as an equation where the final pressure is zero:")
    print("ρ * g * H - (1/2) * ρ * v^2 = 0\n")

    # Step 3: Display the equation with the given numerical values.
    # This is the final equation we need to solve for v.
    print("Substituting the numerical values, the final equation is:")
    print(f"{rho} * {g} * {H} - 0.5 * {rho} * v^2 = 0\n")

    # Step 4: Solve the equation for v.
    # The equation simplifies because ρ can be cancelled from both sides:
    # g * H = 0.5 * v^2
    # Rearranging for v^2: v^2 = 2 * g * H
    v_squared = 2 * g * H
    v = math.sqrt(v_squared)

    # Step 5: Print the step-by-step solution and the final result.
    print("Solving the equation for v:")
    initial_pressure = rho * g * H
    print(f"{initial_pressure} - {0.5 * rho} * v^2 = 0")
    print(f"{initial_pressure} = {0.5 * rho} * v^2")
    print(f"v^2 = {initial_pressure} / {0.5 * rho}")
    print(f"v^2 = {v_squared}")
    print(f"v = sqrt({v_squared})")
    print(f"The required flow speed v is {v:.2f} m/s.")

solve_flow_speed()
print(f"\n<<<14.14>>>")