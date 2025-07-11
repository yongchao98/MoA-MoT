import sympy as sp

def solve_river_crossing_problem():
    """
    Calculates the total downstream drift of a boat with a specific travel pattern.
    """
    # 1. Define the symbolic variables for the problem.
    # L: width of the river
    # v0: maximum flow velocity at the center of the river
    # v: boat's speed relative to the water, perpendicular to the flow
    # y: variable for the distance from the home bank
    L, v0, v, y = sp.symbols('L v_0 v y', positive=True, real=True)

    # 2. Define the river's flow velocity as a function of distance y from the shore.
    # The flow is a parabolic profile: v_flow(y) = k*y*(L-y).
    # We find k by using the condition that v_flow(L/2) = v0.
    # v0 = k * (L/2) * (L - L/2) => v0 = k * L^2 / 4 => k = 4*v0 / L^2
    k = 4 * v0 / L**2
    v_flow = k * y * (L - y)

    # 3. Calculate the drift.
    # The boat's velocity in the y-direction is v. dt = dy / v.
    # The drift rate in the x-direction is dx/dt = v_flow(y).
    # So, the drift per unit change in y is dx/dy = (dx/dt) * (dt/dy) = v_flow / v.
    drift_rate_per_y = v_flow / v

    # 4. Integrate to find the drift for each part of the journey.
    
    # Part 1: Outbound trip from y=0 to y=3L/4.
    turn_around_point_y = 3 * L / 4
    outbound_drift = sp.integrate(drift_rate_per_y, (y, 0, turn_around_point_y))

    # Part 2: Return trip from y=3L/4 to y=0.
    # The drift at any given y is the same, regardless of the direction in the y-axis.
    # The boat travels over the same y-range, so the drift is identical to the outbound trip.
    return_drift = sp.integrate(drift_rate_per_y, (y, 0, turn_around_point_y))
    
    # 5. Calculate the total drift by summing the two parts.
    total_drift = outbound_drift + return_drift

    # 6. Print the results in a step-by-step manner.
    print(f"The river's flow velocity profile is given by: v_flow(y) = {v_flow}")
    print(f"The boat travels from the home bank (y=0) to a turning point at y = 3L/4.")
    print(f"The downstream drift during the outbound trip is: {outbound_drift}")
    print("\nThe boat then returns from y=3L/4 back to the home bank (y=0).")
    print(f"The downstream drift during the return trip is: {return_drift}")
    
    print("\nThe total distance between the starting and returning points is the sum of these two drifts.")
    
    # Extract the numerical coefficients for the final formatted output
    # We divide by the symbolic part to isolate the numerical fraction.
    numerical_part = total_drift / (L * v0 / v)
    num, den = sp.fraction(numerical_part)

    print(f"\nFinal Equation: Distance = {total_drift}")
    print("------------------------------------------")
    print("The numbers in the final equation are:")
    print(f"Numerator coefficient: {num}")
    print(f"Denominator coefficient: {den}")


solve_river_crossing_problem()