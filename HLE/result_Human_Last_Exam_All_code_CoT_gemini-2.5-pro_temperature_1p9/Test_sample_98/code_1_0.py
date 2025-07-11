def solve_tank_shape():
  """
  Determines and prints the shape of the water surface in a half-filled icosahedron tank.
  
  The reasoning is as follows:
  1. An icosahedron has central symmetry.
  2. A tank that is 'half-filled' with water means the water level is at a plane that
     bisects the volume of the tank.
  3. Due to central symmetry, this volume-bisecting plane must pass through the geometric center.
  4. The tank is standing on a face, so the water surface is parallel to this face.
  5. The cross-section of an icosahedron through its center and parallel to a face is a regular hexagon.
  """
  shape = "A regular hexagon"
  print(f"When the icosahedron tank is half-filled, the shape of the water surface will be: {shape}.")

solve_tank_shape()