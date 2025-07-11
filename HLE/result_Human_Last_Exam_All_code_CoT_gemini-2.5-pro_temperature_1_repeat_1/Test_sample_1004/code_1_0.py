import sympy

def solve_river_problem():
    """
    Calculates the total downstream drift of a boat in a river with a triangular flow profile.
    """
    # 1. Define symbolic variables
    # L: River width
    # v0: Maximum flow velocity at the center
    # v: Boat's speed relative to water, perpendicular to the flow
    # y: Distance from the starting bank
    L, v0, v, y = sympy.symbols('L v_0 v y', positive=True, real=True)

    # 2. Define the piecewise flow velocity function
    # The flow is 0 at y=0 and y=L, and maximum (v0) at y=L/2.
    # It increases linearly to the center and then decreases linearly.
    v_flow = sympy.Piecewise(
        (2 * v0 * y / L, y <= L / 2),
        (2 * v0 * (L - y) / L, y > L / 2)
    )

    # 3. Calculate outbound drift (from y=0 to y=3L/4)
    # The time element dt is related to the cross-stream distance dy by dt = dy / v.
    # The drift dx is v_flow * dt. We integrate dx = (v_flow / v) * dy.
    # The integral is over the path from y=0 to y=3L/4.
    outbound_drift = sympy.integrate(v_flow / v, (y, 0, 3 * L / 4))

    # 4. Calculate return drift (from y=3L/4 to y=0)
    # The boat's velocity across the river is now -v, so dt = dy / (-v).
    # We integrate dx = (v_flow / (-v)) * dy over the path from y=3L/4 to y=0.
    return_drift = sympy.integrate(v_flow / (-v), (y, 3 * L / 4, 0))

    # 5. Calculate total drift and simplify
    total_drift = outbound_drift + return_drift
    total_drift_simplified = sympy.simplify(total_drift)

    # 6. Print the results in a clear, step-by-step manner
    print("--- River Boat Drift Calculation ---\n")
    print("Problem Setup:")
    print(f"River width: L")
    print(f"Boat's speed across river: v")
    print(f"Max flow velocity at center: v_0")
    print(f"Flow velocity v_flow(y) = {v_flow}\n")

    print("--- Calculation Steps ---\n")
    print(f"1. Drift during outbound trip (y from 0 to 3L/4):")
    print(f"   delta_x1 = integral(v_flow(y) / v, dy) from 0 to 3L/4")
    print(f"   delta_x1 = {outbound_drift}\n")

    print(f"2. Drift during return trip (y from 3L/4 to 0):")
    print(f"   delta_x2 = integral(v_flow(y) / (-v), dy) from 3L/4 to 0")
    print(f"   delta_x2 = {return_drift}\n")

    print("--- Final Result ---\n")
    print(f"Total downstream drift = Outbound Drift + Return Drift")
    print(f"Total downstream drift = {total_drift}")
    print(f"\nThe simplified final expression for the total distance is:")
    print(f"Distance = {total_drift_simplified}\n")
    
    # 7. Extract the numbers from the final equation as requested
    coeff = total_drift_simplified / (L * v0 / v)
    num, den = sympy.fraction(coeff)
    
    print("--- Final Equation Breakdown ---\n")
    print("The final equation for the total distance is structured as:")
    print(f"Distance = ({num} / {den}) * (L * v_0 / v)")
    print("\nWhere the numbers in the equation are:")
    print(f"Numerator: {num}")
    print(f"Denominator: {den}")


if __name__ == '__main__':
    solve_river_problem()