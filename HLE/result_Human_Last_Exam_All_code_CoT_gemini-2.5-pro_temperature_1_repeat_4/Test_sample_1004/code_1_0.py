import sympy

def solve_boat_problem():
    """
    Calculates the total downstream drift of a boat crossing a river and returning.
    """
    # 1. Define symbolic variables
    # L: width of the river
    # v_0: maximum flow velocity at the center
    # v: boat's speed relative to water, perpendicular to the flow
    # y: distance from the starting bank
    L, v_0, v, y = sympy.symbols('L v_0 v y', positive=True)

    # 2. Define the river's flow velocity u(y)
    # The velocity is a piecewise function. It is proportional to the distance from the
    # nearest shore. The constant of proportionality k is found from the condition
    # that the velocity is v_0 at the center (y=L/2), which gives k = 2*v_0/L.
    flow_velocity_part1 = (2 * v_0 / L) * y
    flow_velocity_part2 = (2 * v_0 / L) * (L - y)
    u = sympy.Piecewise(
        (flow_velocity_part1, y <= L / 2),
        (flow_velocity_part2, y > L / 2)
    )

    # 3. Calculate the downstream drift
    # The boat's velocity components relative to the ground are (dx/dt, dy/dt) = (u(y), v).
    # From this, we can find the drift dx for a small step dy: dt = dy/v.
    # So, dx = u(y) * dt = (u(y)/v) * dy.
    # The integrand for the drift calculation is therefore u(y)/v.
    integrand = u / v

    # The outward journey is from y=0 to y = L - L/4 = 3*L/4.
    drift_outward = sympy.integrate(integrand, (y, 0, 3 * L / 4))

    # The return journey is from y=3L/4 back to y=0. At any given y-position, the
    # river's flow u(y) is the same. Since the boat's speed v across the river is
    # also constant, the time taken to cross any small segment dy is the same for
    # both trips. Therefore, the drift on the return trip is identical to the drift
    # on the outward trip.
    drift_return = drift_outward

    # The total distance from the starting point is the sum of the drifts from both trips.
    total_drift = drift_outward + drift_return

    # 4. Format and print the final answer
    # The result is an expression involving L, v_0, and v. We extract the numerical
    # part to display it clearly as requested. We divide the total_drift by
    # (L * v_0 / v) to get the numerical coefficient as a fraction.
    numerical_coefficient = sympy.simplify(total_drift / (L * v_0 / v))
    num, den = numerical_coefficient.as_numer_denom()

    print("The problem asks for the distance between the boat's returning position on the home bank and its original starting point.")
    print("This distance is the total downstream drift, which is the sum of the drift during the outward and return trips.")
    print("\nThe final equation for the total distance is:")
    # The following print statement shows each number (7 and 8) in the final equation.
    print(f"Distance = ({num} * L * v_0) / ({den} * v)")

if __name__ == '__main__':
    solve_boat_problem()