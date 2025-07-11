import sympy

def solve_boat_problem():
    """
    Calculates the downstream drift of a boat with a specific travel pattern.
    """
    # 1. Define symbolic variables for the problem parameters.
    # L: width of the river
    # v0: maximum flow velocity at the center of the river
    # v: boat's speed relative to the water, perpendicular to the flow
    # y: distance from the starting bank
    L, v0, v, y = sympy.symbols('L v_0 v y', positive=True, real=True)

    # 2. Define the river's flow velocity v_flow(y) as a piecewise function.
    # The velocity is proportional to the distance from the nearest shore.
    # It's 0 at y=0 and y=L, and v0 at y=L/2.
    # For 0 <= y <= L/2, v_flow(y) = k*y. Since v_flow(L/2)=v0, k=2*v0/L.
    # For L/2 < y <= L, v_flow(y) = k*(L-y).
    v_flow_part1 = (2 * v0 / L) * y
    v_flow_part2 = (2 * v0 / L) * (L - y)
    v_flow = sympy.Piecewise(
        (v_flow_part1, y <= L/2),
        (v_flow_part2, y > L/2)
    )

    # 3. Calculate the downstream drift for the outbound trip.
    # The boat travels from y=0 to y=3L/4.
    # The rate of change of drift with respect to y is dx/dy = v_flow(y) / v.
    # We integrate this from 0 to 3L/4.
    integrand = v_flow / v
    x_out = sympy.integrate(integrand, (y, 0, 3 * L / 4))

    # 4. Calculate the total drift.
    # The drift on the return trip (from y=3L/4 to y=0) is identical to the outbound trip.
    # Total drift = drift_out + drift_in = 2 * drift_out.
    total_drift = 2 * x_out

    # 5. Simplify the final expression and print the results clearly.
    final_expression = sympy.simplify(total_drift)

    # Extract the numerical coefficients from the final expression.
    # The expression is of the form (A * ...)/(B * ...)
    num, den = final_expression.as_numer_denom()
    num_coeff = num.as_coeff_Mul()[0]
    den_coeff = den.as_coeff_Mul()[0]

    print("The final equation for the total drift distance is:")
    # Use sympy.pretty_print for a more readable mathematical format
    sympy.pprint(final_expression, use_unicode=True)
    
    print("\nThe numbers in the final equation are:")
    print(f"Numerator: {num_coeff}")
    print(f"Denominator: {den_coeff}")

if __name__ == '__main__':
    solve_boat_problem()