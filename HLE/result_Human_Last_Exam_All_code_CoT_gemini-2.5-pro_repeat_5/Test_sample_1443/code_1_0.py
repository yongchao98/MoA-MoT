import math

def get_asymptote_angles():
    """
    This function calculates and explains the angles of the asymptotes of a specific conic.

    The problem describes a conic passing through five points: A', B', C', O, and H'.
    - A'B'C' is a triangle derived from a base triangle ABC and a point X on its circumcircle.
    - H' is the orthocenter of A'B'C'.
    - O is the circumcenter of ABC.

    Key geometric insights reveal the solution:
    1. A conic through the vertices of a triangle (A'B'C') and its orthocenter (H') is always a rectangular hyperbola.
       This means its asymptotes are perpendicular.
    2. The orientation of this hyperbola depends only on the angle 'delta' between the line 'l' and the line 'BC'. It is
       surprisingly independent of the triangle's specific shape (alpha, beta, gamma) and the choice of point X.

    The final angles of the asymptotes with respect to line BC are found to be -delta and (pi/2 - delta).
    """

    print("The two perpendicular asymptotes make angles with the line BC given by the following equations.")
    print("The angles depend only on 'delta', the angle between line l and line BC.\n")

    # First angle
    # The equation is: theta_1 = -1 * delta
    print("Equation for the first angle, theta_1:")
    print("-1 * delta")
    print("The number in this equation is:")
    print(-1)

    print("-" * 20)

    # Second angle
    # The equation is: theta_2 = (pi / 2) - delta = (1/2) * pi + (-1) * delta
    print("Equation for the second angle, theta_2:")
    print("(1/2) * pi + (-1) * delta")
    print("The numbers in this equation are:")
    print("1, 2, -1")

if __name__ == '__main__':
    get_asymptote_angles()