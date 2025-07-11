import math

def get_asymptote_angles(delta_rad):
    """
    Calculates the angles of the asymptotes with respect to the line BC.

    Args:
        delta_rad: The angle between line BC and line l, in radians.

    Returns:
        A tuple containing the two asymptote angles in radians.
    """
    # Based on the geometric derivation, the angles of the asymptotes of the
    # conic relative to the line BC are 2*delta and 2*delta + pi/2.
    
    angle1_rad = 2 * delta_rad
    angle2_rad = 2 * delta_rad + math.pi / 2
    
    return (angle1_rad, angle2_rad)

def main():
    """
    Main function to explain and display the solution.
    """
    # The problem asks for the angles between the asymptotes of a specific conic
    # and the line BC. The inputs are the triangle angles alpha, beta, gamma,
    # and the angle delta between line BC and an arbitrary line l.
    
    # The solution for the angles of the asymptotes is independent of the
    # triangle's angles (alpha, beta, gamma) and the specific location of
    # point X on the circumcircle. It only depends on the angle delta.
    
    # The two angles, let's call them phi_1 and phi_2, are perpendicular,
    # as the conic is a rectangular hyperbola.
    # The angles with respect to the line BC are given by the expressions:
    
    # phi_1 = 2 * delta
    # phi_2 = 2 * delta + pi/2
    
    print("The angles of the asymptotes with respect to the line BC are given by the expressions:")
    print("Let delta be the angle between line BC and line l.")
    
    # We output the expressions themselves as the answer.
    # The prompt says: "Remember in the final code you still need to output each number
    # in the final equation!". This is interpreted as printing the final expressions.
    
    print("\nIn radians:")
    print("Angle 1 = 2 * delta")
    print("Angle 2 = 2 * delta + pi/2")
    
    print("\nIn degrees:")
    print("Angle 1 = 2 * delta_degrees")
    print("Angle 2 = 2 * delta_degrees + 90")

    # The user is not asked to provide a value for delta,
    # so we will not ask for input and will only print the derived expressions.

if __name__ == "__main__":
    main()
