import sympy
from sympy import symbols, Piecewise, integrate

def solve_boat_drift():
    """
    Calculates the total downstream drift of a boat making a partial round trip across a river.
    """
    # Step 1: Define the symbolic variables for river width (L),
    # max flow velocity (v0), boat's relative speed (v), and position across the river (y).
    L, v_0, v, y = symbols('L v_0 v y', real=True, positive=True)

    # Step 2: Define the piecewise river flow velocity profile.
    # The velocity is proportional to the distance from the nearest shore.
    # It's 0 at y=0 and y=L, and v_0 at the center y=L/2.
    v_flow_part1 = (2 * v_0 / L) * y
    v_flow_part2 = (2 * v_0 / L) * (L - y)
    v_flow = Piecewise(
        (v_flow_part1, y <= L / 2),
        (v_flow_part2, y > L / 2)
    )

    # Step 3: Calculate the drift for the outbound trip.
    # The boat travels from y=0 to y = L - L/4 = 3*L/4.
    # The time to cross a small distance dy is dt = dy / v.
    # The drift for this segment is dx = v_flow(y) * dt = (v_flow(y) / v) * dy.
    # We integrate this from y=0 to y=3*L/4.
    drift_outbound = integrate(v_flow / v, (y, 0, 3 * L / 4))

    # Step 4: Calculate the total drift.
    # The return trip from y=3*L/4 to y=0 covers the same y-distances
    # and thus experiences the same flow profile, resulting in an identical drift.
    total_drift = 2 * drift_outbound

    # Step 5: Display the final result and the numbers in the equation.
    print("The final equation for the total downstream distance is:")

    # Extract the numerical coefficient from the symbolic result.
    # The symbolic part of the expression is (L * v_0 / v)
    symbolic_part = (L * v_0) / v
    numeric_coefficient = total_drift / symbolic_part
    
    # Get the numerator and denominator of the coefficient
    num = sympy.numer(numeric_coefficient)
    den = sympy.denom(numeric_coefficient)

    print(f"Distance = ({num} * L * v_0) / ({den} * v)")
    print("\nBreaking down the equation:")
    print(f"The number in the numerator is: {num}")
    print(f"The number in the denominator is: {den}")
    print(f"The other variables are L (river width), v_0 (max flow speed), and v (boat's speed).")


# Execute the function to solve the problem
solve_boat_drift()