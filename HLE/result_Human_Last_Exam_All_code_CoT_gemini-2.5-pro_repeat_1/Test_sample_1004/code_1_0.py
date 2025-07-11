import sympy

def solve_river_crossing_problem():
    """
    This function calculates the total downstream displacement of a boat
    crossing a river with a parabolic velocity profile.
    """
    # Step 1: Define symbolic variables for the problem.
    # L: width of the river
    # v_0: maximum flow velocity at the center of the river
    # v: the boat's speed relative to the water, perpendicular to the flow
    # y: the distance from the starting bank
    L, v_0, v, y = sympy.symbols('L v_0 v y', real=True, positive=True)

    # Step 2: Formulate the river's flow velocity profile u(y).
    # The profile is parabolic: u(y) = k*y*(L-y).
    # We know u(L/2) = v_0, so we can find k.
    # v_0 = k * (L/2) * (L - L/2) => k = 4*v_0/L**2
    u_y = (4 * v_0 / L**2) * y * (L - y)

    # Step 3: Calculate the downstream drift for the trip out.
    # The boat travels from y=0 to y=3L/4.
    # The boat's velocity across the river is dy/dt = v, so dt = dy/v.
    # The drift dx = u(y) * dt = (u(y)/v) * dy.
    # We integrate to find the total drift for this leg.
    integrand_out = u_y / v
    x_out = sympy.integrate(integrand_out, (y, 0, 3 * L / 4))

    # Step 4: Calculate the downstream drift for the trip back.
    # The boat travels from y=3L/4 to y=0.
    # The boat's velocity across the river is dy/dt = -v, so dt = -dy/v.
    # The drift dx = u(y) * dt = -(u(y)/v) * dy.
    # We integrate from y=3L/4 to y=0.
    integrand_back = -u_y / v
    x_back = sympy.integrate(integrand_back, (y, 3 * L / 4, 0))

    # Step 5: Sum the drifts to find the total distance.
    total_distance = sympy.simplify(x_out + x_back)
    
    # Extract the numerical parts of the final expression to print them clearly.
    # The expression is of the form: (numerator/denominator) * (v_0 * L / v)
    coefficient = sympy.simplify(total_distance / (v_0 * L / v))
    num, den = coefficient.as_numer_denom()

    print("The final equation for the total downstream distance is:")
    print(f"Distance = ({num} * v_0 * L) / ({den} * v)")
    print("\nIn this equation:")
    print(f"- The numerator of the coefficient is: {num}")
    print(f"- The denominator of the coefficient is: {den}")
    print("- v_0 is the maximum flow velocity.")
    print("- L is the width of the river.")
    print("- v is the boat's speed relative to the water.")

if __name__ == '__main__':
    solve_river_crossing_problem()