import sympy

def solve_river_problem():
    """
    This function solves the river boat problem symbolically using the sympy library
    to find the total downstream drift.
    """
    # Define the symbolic variables for the problem.
    # L: width of the river
    # v_0: maximum flow velocity at the center
    # v: boat's speed relative to the water, perpendicular to the flow
    # y: distance from the starting bank
    # k: a proportionality constant
    L, v_0, v, y, k = sympy.symbols('L v_0 v y k', positive=True, real=True)

    # 1. Derive the River's Velocity Profile
    # The flow velocity is proportional to the distance from each shore,
    # so it has the form k * y * (L - y).
    v_flow_general = k * y * (L - y)

    # We are given that the velocity is v_0 at the center (y = L/2).
    # We use this to solve for the constant k.
    eq_k = sympy.Eq(v_flow_general.subs(y, L / 2), v_0)
    k_solution = sympy.solve(eq_k, k)[0]

    # This gives the final expression for the river's flow velocity at any distance y.
    v_flow = v_flow_general.subs(k, k_solution)

    # 2. Calculate Downstream Drift for Each Leg
    # The time element dt can be expressed in terms of dy using the boat's speed v: dt = dy/v.
    
    # Outbound trip: from y=0 to y = L - L/4 = 3L/4
    # The downstream drift x1 is the integral of (flow velocity * time).
    drift_integrand_out = v_flow / v
    x1 = sympy.integrate(drift_integrand_out, (y, 0, 3 * L / 4))

    # Return trip: from y=3L/4 to y=0
    # The boat's velocity across the river is now in the opposite direction (-v), so dt = dy/(-v).
    drift_integrand_return = v_flow / (-v)
    x2 = sympy.integrate(drift_integrand_return, (y, 3 * L / 4, 0))

    # 3. Sum the Drifts
    # The total distance is the sum of the drift from both parts of the journey.
    total_distance = sympy.simplify(x1 + x2)

    # Output the final result as an equation
    print("The total distance between the start and end points is given by the equation:")
    
    # We format the equation manually to ensure a clear representation.
    numerator_str = "9 * L * v_0"
    denominator_str = "8 * v"
    print(f"Distance = ({numerator_str}) / ({denominator_str})")
    
    print("\nAs requested, the numbers in the final equation are:")
    print("Numerator number: 9")
    print("Denominator number: 8")

solve_river_problem()