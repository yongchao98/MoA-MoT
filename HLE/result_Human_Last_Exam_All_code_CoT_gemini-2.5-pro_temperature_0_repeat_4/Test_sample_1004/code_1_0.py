import sympy

def solve_river_crossing_problem():
    """
    This script symbolically solves the river crossing problem to find the total
    downstream drift of the boat.
    """
    # Step 1: Define the symbolic variables for the problem.
    # L: width of the river
    # v_0: maximum flow velocity at the center of the river
    # v: boat's speed relative to the water, perpendicular to the flow
    # y: variable for distance from the starting bank
    # k: a proportionality constant for the flow profile
    L, v_0, v, y, k = sympy.symbols('L v_0 v y k')

    print("Step 1: Defining the river's velocity profile u(y).")
    # The flow velocity is proportional to the distance from both shores.
    # This can be modeled with a parabolic equation: u(y) = k * y * (L - y).
    u_expr = k * y * (L - y)

    # We find the constant 'k' using the condition that velocity is v_0 at the center (y = L/2).
    # We solve the equation: v_0 = k * (L/2) * (L - L/2) for k.
    k_value = sympy.solve(u_expr.subs(y, L/2) - v_0, k)[0]

    # Substitute the value of k back into the velocity profile expression.
    u_y = u_expr.subs(k, k_value)
    print(f"The derived velocity profile of the river is: u(y) = {sympy.simplify(u_y)}\n")

    # The boat's velocity perpendicular to the flow is v. The time element dt = dy / v.
    # The downstream drift dx is u(y) * dt.
    # To find the total drift for a segment, we integrate u(y)/v with respect to y.
    integrand = u_y / v

    print("Step 2: Calculating the downstream drift for the outbound journey.")
    # The boat travels from y=0 to y = L - L/4 = 3L/4.
    # We integrate from y=0 to y=3L/4 to find the drift x_1.
    x_1 = sympy.integrate(integrand, (y, 0, 3*L/4))
    print(f"The drift during the outbound trip to y=3L/4 is: {sympy.simplify(x_1)}\n")

    print("Step 3: Calculating the downstream drift for the inbound journey.")
    # The boat travels back from y=3L/4 to y=0.
    # The velocity component is -v, so dt = dy/(-v). The integration is from 3L/4 to 0.
    # The integral is ∫[3L/4 to 0] u(y)/(-v) dy, which simplifies to ∫[0 to 3L/4] u(y)/v dy.
    # This is the same as the outbound drift.
    x_2 = sympy.integrate(u_y / (-v), (y, 3*L/4, 0))
    print(f"The drift during the inbound trip from y=3L/4 is: {sympy.simplify(x_2)}\n")

    print("Step 4: Calculating the total downstream drift.")
    # The total distance is the sum of the two drifts.
    total_drift = x_1 + x_2
    simplified_drift = sympy.simplify(total_drift)
    print(f"Total Drift = (Drift Out) + (Drift In) = {sympy.simplify(x_1)} + {sympy.simplify(x_2)}")
    print(f"The final simplified expression for the total drift is: {simplified_drift}\n")

    print("--- Final Answer ---")
    print("The distance between the boat's returning position and its original starting point is given by the equation:")
    # Extract the numerical coefficients for the final output as requested.
    # The result is of the form (Num/Den) * v_0 * L / v
    num, den = simplified_drift.as_numer_denom()
    # We assume the symbolic part is v_0*L/v
    coeff = simplified_drift / (v_0 * L / v)
    
    print(f"Total Drift = ({coeff}) * (v_0 * L / v)")
    print("\nThe numbers in the final equation are:")
    print(f"Numerator coefficient: {coeff.p}")
    print(f"Denominator coefficient: {coeff.q}")

if __name__ == '__main__':
    solve_river_crossing_problem()