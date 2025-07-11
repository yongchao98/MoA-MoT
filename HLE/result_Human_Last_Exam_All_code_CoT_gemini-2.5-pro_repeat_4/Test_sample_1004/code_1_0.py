import sympy

def solve_river_problem():
    """
    Calculates the total downstream drift of a boat crossing a river and returning.
    """
    # Define the symbols used in the problem. 'positive=True' helps sympy with simplifications.
    L, v0, v, y = sympy.symbols('L v_0 v y', positive=True, real=True)

    print("Step 1: Define the river flow velocity u(y).")
    print("The flow is modeled as a piecewise linear function.")
    # For the first half of the river (0 <= y <= L/2)
    u_first_half = (2 * v0 / L) * y
    # For the second half of the river (L/2 <= y <= L)
    u_second_half = (2 * v0 / L) * (L - y)
    
    print(f"  - For 0 <= y <= L/2, the flow velocity is u(y) = {sympy.pretty(u_first_half, use_unicode=False)}")
    print(f"  - For L/2 <= y <= L, the flow velocity is u(y) = {sympy.pretty(u_second_half, use_unicode=False)}\n")

    # The boat's speed perpendicular to the flow is v.
    # The time to cross a small distance dy is dt = dy / v.
    # The drift dx during this time is dx = u(y) * dt = u(y) * (dy / v).

    print("Step 2: Calculate the downstream drift for the outbound trip (from y=0 to y=3L/4).")
    # This trip crosses the center, so we must split the integral into two parts.
    # Part 1.1: Drift from y=0 to y=L/2
    drift_out_1 = sympy.integrate(u_first_half / v, (y, 0, L / 2))
    # Part 1.2: Drift from y=L/2 to y=3L/4
    drift_out_2 = sympy.integrate(u_second_half / v, (y, L / 2, 3 * L / 4))
    
    # Total drift for the outbound trip
    total_drift_out = drift_out_1 + drift_out_2
    print(f"The drift during the outbound trip is the sum of two integrals, resulting in: {sympy.pretty(total_drift_out, use_unicode=False)}\n")

    print("Step 3: Calculate the downstream drift for the return trip (from y=3L/4 to y=0).")
    # The magnitude of the drift during the return trip is identical to the outbound trip
    # because the boat spends the same amount of time at each y-position.
    total_drift_return = total_drift_out
    print(f"The drift during the return trip is the same: {sympy.pretty(total_drift_return, use_unicode=False)}\n")

    print("Step 4: Calculate the total distance from the starting point.")
    # The total drift is the sum of the drift from the outbound and return trips.
    total_drift = total_drift_out + total_drift_return
    print("The total distance is the sum of the drifts from both trips.")
    print(f"Total Drift = {sympy.pretty(total_drift_out, use_unicode=False)} + {sympy.pretty(total_drift_return, use_unicode=False)}")
    
    # Extract the coefficients for the final print statement
    num, den = sympy.fraction(total_drift)
    num_coeff = num.as_coeff_mul(L, v0)[0]
    den_coeff = den.as_coeff_mul(v)[0]
    
    print("\n--- Final Answer ---")
    print(f"The final equation for the total drift distance is: {sympy.pretty(total_drift, use_unicode=False)}")
    print(f"This can be written as: ({num_coeff} * L * v_0) / ({den_coeff} * v)")
    print(f"The numbers in the final equation are the numerator {num_coeff} and the denominator {den_coeff}.")

if __name__ == '__main__':
    solve_river_problem()