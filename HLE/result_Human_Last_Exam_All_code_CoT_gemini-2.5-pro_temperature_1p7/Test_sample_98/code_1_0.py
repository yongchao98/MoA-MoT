def solve_icosahedron_problem():
    """
    Explains the reasoning to determine the shape of the water surface
    in a half-filled icosahedron tank standing on one face.
    """
    
    print("Here is the step-by-step reasoning to find the shape of the water surface:")
    
    print("\nStep 1: Analyze the geometry of the tank.")
    print("The tank is an icosahedron, a regular polyhedron with 20 identical triangular faces and 12 vertices.")
    print("When it stands on one face, that face is horizontal. Due to the icosahedron's symmetry, there is a parallel face at the top.")
    
    print("\nStep 2: Determine the water level.")
    print("The tank is 'half-filled', meaning the water occupies half of the total volume.")
    print("Because of the symmetry between the top and bottom faces, the plane that divides the volume in half is the horizontal plane passing through the geometric center of the icosahedron.")
    
    print("\nStep 3: Identify the shape of the water's surface.")
    print("The water's surface corresponds to the shape of the cross-section created by this central plane.")
    print("Of the 12 vertices of the icosahedron, 3 form the base and 3 form the top.")
    print("The other 6 vertices lie exactly on the central plane, halfway between the top and bottom.")
    
    print("\nConclusion:")
    print("These 6 vertices on the central plane form a regular polygon. Therefore, the shape of the water surface is a regular hexagon.")

solve_icosahedron_problem()
