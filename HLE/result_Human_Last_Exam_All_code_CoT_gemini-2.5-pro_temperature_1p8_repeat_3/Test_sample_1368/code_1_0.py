def solve_inscribed_square():
    """
    This function determines the shortest sequence of commands to construct a square
    inscribed in a circle.

    The plan is as follows:
    Let the center of the circle be O and a point on the circumference be A.

    1.  (L)ine: Draw a line through O and A. This line is a diameter and
        intersects the circle at a new point, C. We now have one diagonal, AC.

    2.  (C)ircle: Draw a circle centered at A with radius AC.

    3.  (C)ircle: Draw a circle centered at C with radius CA. These two circles
        intersect at two new points, P and Q.

    4.  (L)ine: Draw a line through P and Q. This line is the perpendicular
        bisector of the diameter AC and forms the second diameter, BD.

    The vertices of the square are A, B, C, and D.
    The sequence of commands is LCCL.
    """
    sequence = "LCCL"
    print(sequence)

solve_inscribed_square()