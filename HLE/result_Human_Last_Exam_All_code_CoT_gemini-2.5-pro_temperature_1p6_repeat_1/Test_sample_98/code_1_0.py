def get_water_surface_shape():
  """
  Determines the shape of the water surface in a half-filled icosahedron tank
  standing on one of its faces.

  Reasoning:
  1. An icosahedron is a regular polyhedron with 20 equilateral triangular faces.
  2. When standing on a face, the base is a horizontal triangle. Due to symmetry, there is a parallel triangular face at the top.
  3. A half-filled tank means the water occupies half the total volume.
  4. The plane that bisects the volume of the icosahedron is the horizontal plane passing through its geometric center, which is located at mid-height.
  5. The shape of the water's surface is therefore the shape of the cross-section of the icosahedron at its center, parallel to the base face.
  6. This cross-section is a regular hexagon. This is because the cross-section must respect the 3-fold rotational symmetry (around the axis connecting the top and bottom faces) and the central inversion symmetry of the icosahedron.
  """
  shape = "A regular hexagon"
  print(f"When the icosahedron tank is half-filled, the shape of the water surface will be: {shape}")

get_water_surface_shape()