def solve_boat_drift():
    """
    Calculates and prints the formula for the boat's total downstream drift.

    The problem involves two parts of a journey:
    1. Outward trip: from y=0 to y=3L/4
    2. Return trip: from y=3L/4 to y=0

    The river's velocity profile is u(y) = (4*v_0/L^2) * y * (L-y).
    The boat's speed across the river is a constant, v.

    The drift for a small dy is dx = u(y) * dt, where dt = dy/v.
    So, dx = u(y)/v * dy.

    The drift for the first leg is the integral of dx from 0 to 3L/4:
    Drift_1 = integral( (4*v_0/(v*L^2)) * (L*y - y^2) dy ) from 0 to 3L/4
            = (4*v_0/(v*L^2)) * [L*y^2/2 - y^3/3] from 0 to 3L/4
            = (4*v_0/(v*L^2)) * (L*(3L/4)^2/2 - (3L/4)^3/3)
            = (4*v_0/(v*L^2)) * (9*L^3/32 - 9*L^3/64)
            = (4*v_0/(v*L^2)) * (9*L^3/64)
            = (9 * v_0 * L) / (16 * v)

    The drift for the return leg has the same magnitude because the boat covers the
    same y-span, and the river flow u(y) depends only on y.
    Drift_2 = (9 * v_0 * L) / (16 * v)

    The total drift is the sum of the two parts.
    Total Drift = Drift_1 + Drift_2
    """

    # Symbolic representation of variables
    L = "L"   # Width of the river
    v0 = "v_0" # Maximum flow velocity at the center
    v = "v"   # Boat's speed relative to the water

    # The final formula is derived from adding the drift from both journey parts.
    # Total Drift = 2 * (9 * v_0 * L) / (16 * v)
    # This simplifies to (18 / 16) * (v_0 * L / v)
    
    # Coefficients of the final simplified formula
    numerator_coefficient = 9
    denominator_coefficient = 8

    # Print the final formula for the total distance
    print("The final formula for the distance between the boat's returning position and its original starting point is:")
    print(f"Distance = ({numerator_coefficient} * {v0} * {L}) / ({denominator_coefficient} * {v})")
    
solve_boat_drift()