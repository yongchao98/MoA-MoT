import sympy

def solve_river_crossing_problem():
    """
    Calculates the downstream distance a boat travels on a round trip,
    turning around at a specific point.
    """
    # Step 1: Define symbolic variables
    # L: width of the river
    # v_0: maximum flow velocity at the river's center
    # v: boat's speed relative to the water, perpendicular to the flow
    # y: variable for the position across the river
    L, v_0, v, y, k = sympy.symbols('L v_0 v y k')

    print("### Step 1: Determine the River's Flow Velocity Profile ###")
    print("Let the river width be L. The y-axis points across the river from y=0 to y=L.")
    print("The flow velocity is 0 at y=0 and y=L, and max (v_0) at y=L/2.")
    print("This can be modeled with a parabolic function: v_flow(y) = k * y * (L - y)")

    # Define the base velocity expression
    v_flow_expr = k * y * (L - y)

    # Solve for the constant k using the condition v_flow(L/2) = v_0
    k_eq = sympy.Eq(v_flow_expr.subs(y, L/2), v_0)
    k_sol = sympy.solve(k_eq, k)[0]

    # Final expression for the flow velocity
    v_flow = v_flow_expr.subs(k, k_sol)
    print(f"\nBy solving for k, the velocity profile is: v_flow(y) = {sympy.simplify(v_flow)}\n")

    # Step 2: Calculate downstream drift on the outward journey
    print("### Step 2: Analyze the Outward Journey ###")
    print("The boat travels from y=0 to y = L - L/4 = 3L/4.")
    print("The time to travel across the river is dt = dy / v.")
    print("The downstream drift (dx) is v_flow(y) * dt = (v_flow(y) / v) * dy.")
    print("Integrating this from y=0 to y=3L/4 gives the drift for the first leg (x1).")

    integrand = v_flow / v
    x1 = sympy.integrate(integrand, (y, 0, 3*L/4))
    print(f"\nThe calculated drift for the outward journey is: x1 = {sympy.simplify(x1)}\n")

    # Step 3: Calculate downstream drift on the return journey
    print("### Step 3: Analyze the Return Journey ###")
    print("The boat turns around and travels from y=3L/4 back to y=0.")
    print("The speed across the river is now in the negative y direction, but the magnitude is the same.")
    print("The downstream drift calculation results in the same integral.")

    x2 = x1  # The integral is identical
    print(f"The calculated drift for the return journey is also: x2 = {sympy.simplify(x2)}\n")

    # Step 4: Calculate the total distance
    print("### Step 4: Calculate the Total Distance ###")
    print("The total distance from the start is the sum of the drifts from both legs.")
    total_distance = x1 + x2
    simplified_distance = sympy.simplify(total_distance)
    print(f"Total Distance = x1 + x2 = {simplified_distance}\n")
    
    # Final result breakdown as requested
    print("--- Final Equation ---")
    num, den = sympy.fraction(simplified_distance / (L*v_0/v))
    print(f"Distance = ({num} * L * v_0) / ({den} * v)")
    print(f"The numbers in the final equation are:")
    print(f"  - Numerator coefficient: {num}")
    print(f"  - Denominator coefficient: {den}")


if __name__ == '__main__':
    solve_river_crossing_problem()