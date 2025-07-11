def solve_tank_problem():
    """
    This function determines and prints the shape of the water surface
    in a half-filled icosahedron tank standing on one face.
    """
    
    # An icosahedron is a centrally symmetric polyhedron.
    # When it's half-full, the water level is a plane passing through its center.
    # When an icosahedron stands on one face, this central plane is parallel to the base.
    # The cross-section of an icosahedron through its center and parallel to a face is a regular hexagon.
    
    shape = "A regular hexagon"
    
    print(f"The shape of the water surface will be: {shape}")

solve_tank_problem()