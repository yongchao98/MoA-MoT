import math

def solve_flow_speed():
    """
    Calculates the flow speed 'v' at which the pressure at the bottom of a river becomes zero.
    """
    # Given constants
    rho = 1000  # Density of water in kg/m^3
    g = 10      # Acceleration due to gravity in m/s^2
    H = 10      # Depth of the river in meters

    # --- Explanation ---

    # Step 1: Calculate the initial static gauge pressure at the bottom of the river.
    # This pressure is due to the weight of the water column.
    # P_static_gauge = ρ * g * H
    p_static_gauge = rho * g * H
    print("The initial gauge pressure at the bottom is given by P_static = ρ * g * H.")
    print(f"P_static = {rho} * {g} * {H} = {p_static_gauge} N/m^2.\n")

    # Step 2: Formulate the pressure equation for flowing water.
    # When water flows, the static pressure is reduced by the dynamic pressure term (1/2 * ρ * v^2).
    # The final gauge pressure P_final = P_static - P_dynamic.
    print("When the water flows, the final gauge pressure (P_final) is:")
    print("P_final = (ρ * g * H) - (1/2 * ρ * v^2)\n")

    # Step 3: Set the condition for zero pressure.
    # We want to find the speed v for which P_final = 0.
    print("Setting the final pressure to 0 to find the required speed v:")
    print("0 = (ρ * g * H) - (1/2 * ρ * v^2)\n")

    # Step 4: Solve the equation for v.
    # 1/2 * ρ * v^2 = ρ * g * H
    # v^2 = 2 * g * H
    # v = sqrt(2 * g * H)
    print("Rearranging the equation to solve for v:")
    print("1/2 * ρ * v^2 = ρ * g * H")
    print("v = sqrt(2 * g * H)\n")

    # Step 5: Substitute the values and calculate the result.
    v_squared = 2 * g * H
    v = math.sqrt(v_squared)

    print("Substituting the numerical values into the final equation:")
    print(f"v = sqrt(2 * {g} * {H})")
    print(f"v = sqrt({v_squared})")
    print(f"The required flow speed is {v:.2f} m/s.")

    # The final answer in the required format
    print(f"<<<{v:.2f}>>>")

# Run the solver
solve_flow_speed()