import math

def solve_rotational_symmetry():
    """
    This function determines the rotational symmetry of the provided tiling image.
    
    Step-by-step reasoning:
    1.  The tiling consists of several repeating shapes, notably blue hexagons and tan squares.
        These regular-looking polygons are candidates for centers of rotational symmetry.

    2.  Hypothesis 1: 6-fold symmetry based on the blue hexagons.
        A full circle is 360 degrees. For 6-fold symmetry, the minimum rotation angle that
        leaves the pattern unchanged would be 360 / 6 = 60 degrees.
        Visually inspecting the pattern, a 60-degree rotation around the center of any
        blue hexagon maps the hexagon and all its surrounding tiles perfectly onto
        the positions of identical tiles. This symmetry appears to hold true for the entire pattern.

    3.  Hypothesis 2: 4-fold symmetry based on the tan squares.
        For 4-fold symmetry, the rotation angle would be 360 / 4 = 90 degrees.
        Let's test this. A 90-degree rotation around the center of a square maps the square
        onto itself. However, looking closely at the tiles meeting at the square's corners reveals
        that the arrangement is not 4-fold symmetric. Each corner of a square adjoins one
        yellow pentagon and one orange rhombus. A 90-degree rotation would require a
        pentagon to be mapped onto a rhombus, which is not the case.
        Therefore, the tiling does not have 4-fold rotational symmetry.

    4.  Conclusion: The highest order of rotational symmetry present in the tiling is 6-fold.
    """
    
    # The order of rotational symmetry is the number of times the tiling maps
    # onto itself during a 360-degree rotation.
    order_of_symmetry = 6
    
    # A full rotation is 360 degrees.
    full_circle_degrees = 360
    
    # Calculate the minimum angle of rotation for this symmetry.
    min_rotation_angle = full_circle_degrees / order_of_symmetry
    
    print("The rotational symmetry of the tiling is 6-fold.")
    print("This conclusion is reached by observing that the pattern repeats every 60 degrees when rotated around the center of a blue hexagon.")
    print("\nThe mathematical relationship between the order of symmetry (n), and the minimum angle of rotation is:")
    print("Angle = 360 / n")
    print("\nFor this tiling, the equation with its numbers is:")
    # Printing each number in the final equation
    print(f"{full_circle_degrees} / {order_of_symmetry} = {int(min_rotation_angle)}")

solve_rotational_symmetry()