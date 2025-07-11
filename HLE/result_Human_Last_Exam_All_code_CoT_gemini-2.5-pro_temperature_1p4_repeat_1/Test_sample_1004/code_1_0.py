import sympy as sp

def solve_river_problem():
    """
    Calculates the total downstream drift of a boat crossing a river and returning.
    """
    # 1. Define the symbolic variables
    L, v_0, v, x = sp.symbols('L v_0 v x', real=True, positive=True)

    print("Step 1: Define the river's flow velocity v_flow(x).")
    # 2. Define the piecewise flow velocity profile
    # For 0 <= x <= L/2
    v_flow_part1 = (2 * v_0 / L) * x
    # For L/2 < x <= L
    v_flow_part2 = (2 * v_0 / L) * (L - x)
    print(f"  - For 0 <= x <= L/2: v_flow = {v_flow_part1}")
    print(f"  - For L/2 < x <= L:  v_flow = {v_flow_part2}\n")

    # The boat turns around at x = L - L/4 = 3L/4
    turn_point = 3 * L / 4

    print(f"Step 2: Calculate the downstream drift for a one-way trip to x = {sp.pretty(turn_point)}.")
    # 3. Calculate the integral of the flow velocity for one trip.
    # The drift is Integral(v_flow(x) dt) = Integral(v_flow(x) * dx/v).
    # We need to split the integral at x = L/2.

    # Integral from x=0 to x=L/2
    integral_1 = sp.integrate(v_flow_part1, (x, 0, L/2))
    print(f"  - Drift integral over [0, L/2]: {integral_1}")

    # Integral from x=L/2 to x=3L/4
    integral_2 = sp.integrate(v_flow_part2, (x, L/2, turn_point))
    print(f"  - Drift integral over [L/2, 3L/4]: {integral_2}")

    # Total drift for the outbound trip
    total_integral_one_way = integral_1 + integral_2
    drift_out = total_integral_one_way / v
    print(f"  - Total drift for the outbound trip is: ({total_integral_one_way}) / v = {drift_out}\n")

    # 4. Calculate the total drift for the round trip.
    # The drift on the return journey is identical to the outbound journey.
    print("Step 3: Calculate the total downstream distance for the round trip.")
    total_drift = 2 * drift_out
    
    # Simplify the final expression
    final_distance = sp.simplify(total_drift)
    print("  - The drift on the return trip is the same, so we multiply by 2.")
    print(f"  - Total distance = 2 * ({drift_out}) = {total_drift}\n")

    print("--- Final Result ---")
    print(f"The simplified expression for the total distance is: {final_distance}")

    # 5. Extract and display the numerical coefficients of the final formula
    print("\nTo show the numbers in the final equation:")
    # Get numerator and denominator of the expression
    numerator, denominator = final_distance.as_numer_denom()
    # Get the numerical coefficient from the numerator and denominator
    num_coeff = numerator.as_coeff_Mul()[0]
    den_coeff = denominator.as_coeff_Mul()[0]

    print(f"The final distance can be written as: ({num_coeff} * v_0 * L) / ({den_coeff} * v)")


if __name__ == "__main__":
    solve_river_problem()
