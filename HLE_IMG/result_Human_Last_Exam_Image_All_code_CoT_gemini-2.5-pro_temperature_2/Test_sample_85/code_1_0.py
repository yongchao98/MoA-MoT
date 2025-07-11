import math

def solve_cone_distance():
  """
  Calculates and explains the furthest distance from a point P on the cone's surface.

  The problem states:
  - Base diameter = d
  - Slant height (distance from apex to P) = d

  The solution involves unrolling the cone and finding the furthest point on the flattened 2D surface.
  """
  
  # Step 1: Define cone properties in terms of a symbolic variable 'd'
  # d_val represents a numerical example for d, but the final answer is symbolic.
  # Let's use d=10 for any illustrative calculations, but the logic holds for any d.
  d_val = 10
  
  base_diameter = 'd'
  base_radius = 'd / 2'
  slant_height = 'd'

  # Step 2: Unroll the cone.
  # The unrolled surface is a sector of a circle.
  # Radius of the sector is the slant height.
  sector_radius = slant_height # 'd'
  
  # Arc length of the sector is the base circumference of the cone.
  # Base circumference C = 2 * pi * r = 2 * pi * (d/2) = pi * d
  arc_length = 'pi * d'
  
  # Step 3: Calculate the angle of the sector.
  # Angle (theta) = Arc Length / Radius
  # theta = (pi * d) / d = pi
  theta_radians = math.pi
  theta_degrees = 180
  
  print("Step 1: The cone's slant height L is d, and its base radius r is d/2.")
  print("Step 2: Unrolling the cone's surface gives a sector of a circle with radius L = d.")
  print(f"Step 3: The sector's angle is (base circumference) / (slant height) = (π*d) / d = π radians ({theta_degrees} degrees).")
  print("This means the unrolled surface is a semicircle of radius d.")
  
  # Step 4: Find the furthest point.
  # To find the furthest geodesic distance, we map the cone surface onto a full circle of radius d (two semicircles joined).
  # We place our point P on the edge of this circle, e.g., at (d, 0).
  # The furthest point Q in the circle is the diametrically opposite point, (-d, 0).
  # The distance is the diameter of this circle.
  
  print("Step 4: The furthest point on the surface is found by considering the path from P, through the apex, to the diametrically opposite point on the base.")
  
  # Step 5: The calculation
  # The distance is the path up to the apex (length d) plus the path down to the opposite side (length d).
  # Total distance = d + d = 2d.
  
  coefficient = 2
  
  print("\nThe final equation for the furthest distance is:")
  print(f"Distance = {coefficient} * d")

solve_cone_distance()