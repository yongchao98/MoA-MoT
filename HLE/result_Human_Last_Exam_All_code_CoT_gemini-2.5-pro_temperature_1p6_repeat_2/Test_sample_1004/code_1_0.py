import sympy

def solve_river_crossing_problem():
    """
    This function solves the river crossing problem symbolically using sympy.
    """
    # 1. Define symbolic variables for the problem.
    # L: width of the river
    # v0: maximum flow velocity at the center
    # v: boat's speed relative to water, perpendicular to the flow
    # y: distance from the home bank
    # All are assumed to be positive real numbers.
    L, v0, v, y = sympy.symbols('L v_0 v y', real=True, positive=True)

    # 2. Define the river flow velocity profile u(y).
    # The profile is parabolic: u(y) = k*y*(L-y).
    # We find the constant k from the condition that the velocity is v0 at y=L/2.
    # v0 = k*(L/2)*(L - L/2) = k*L**2/4  => k = 4*v0/L**2
    k = 4 * v0 / L**2
    u_y = k * y * (L - y)

    # 3. Calculate the downstream drift.
    # The boat's velocity component across the river is v_y = v (constant).
    # The time element is dt = dy / v.
    # The drift element is dx = u(y) * dt = (u(y)/v) * dy.

    # Part 1: Outward journey from y=0 to y=3L/4
    # We integrate dx from y=0 to y=3L/4 to find the drift x1.
    integrand = u_y / v
    y_turnaround = 3 * L / 4
    x1 = sympy.integrate(integrand, (y, 0, y_turnaround))

    # Part 2: Return journey from y=3L/4 to y=0
    # The magnitude of the velocity component across the river is still v.
    # The time taken to cover the same distance 3L/4 is the same.
    # The river velocity profile traversed is the same.
    # Therefore, the drift on the return journey, x2, is equal to x1.
    x2 = x1

    # 4. Calculate the total downstream distance from the start point.
    total_distance = x1 + x2

    # 5. Print the results clearly.
    print("The final equation for the total downstream distance is:")
    # Using sympy.pretty_print for a clean mathematical representation.
    sympy.pretty_print(sympy.Eq(sympy.Symbol('Distance'), total_distance))
    
    # Per the instruction, output each number in the final equation.
    # The result is of the form (numerator * v0 * L) / (denominator * v).
    # We can extract the numerical coefficients from the symbolic expression.
    coeff = total_distance / (v0 * L / v)
    
    # sympy.fraction extracts the numerator and denominator of a rational expression.
    numerator_coeff, denominator_coeff = sympy.fraction(coeff)

    print("\nIn the final equation: Distance = (N * v_0 * L) / (D * v)")
    print(f"The number in the numerator is N = {numerator_coeff}")
    print(f"The number in the denominator is D = {denominator_coeff}")

if __name__ == '__main__':
    solve_river_crossing_problem()