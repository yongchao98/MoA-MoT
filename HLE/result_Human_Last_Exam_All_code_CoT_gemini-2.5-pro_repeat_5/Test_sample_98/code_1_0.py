def solve_icosahedron_problem():
    """
    This function explains and determines the shape of the water surface
    in a half-filled icosahedron tank standing on one of its faces.
    """
    print("Here is the step-by-step reasoning to determine the shape of the water surface:")
    print("-" * 70)

    # Step 1: Analyze the geometry of the icosahedron
    print("Step 1: Analyze the geometry of the icosahedron.")
    print("An icosahedron is a polyhedron with 20 identical equilateral triangular faces.")
    print("Crucially, it possesses a center of symmetry. This means that for any point on its surface, the point directly opposite through its geometric center is also on the surface.")
    print("\n")

    # Step 2: Interpret the setup
    print("Step 2: Understand the problem's setup.")
    print("The tank is 'standing on one of its faces', which means this face is the horizontal base. Due to the icosahedron's symmetry, the face directly opposite the base is also horizontal, acting as the top.")
    print("The tank is 'half-filled', meaning the water occupies exactly 50% of the tank's total volume.")
    print("\n")

    # Step 3: Locate the water level
    print("Step 3: Locate the water level.")
    print("Because the icosahedron is centrally symmetric, any plane passing through its center will divide its volume into two equal halves.")
    print("The water's surface is a horizontal plane. For the volume to be exactly half, this plane must pass through the icosahedron's center of symmetry.")
    print("\n")

    # Step 4: Determine the shape of the water surface
    print("Step 4: Determine the shape of the resulting water surface.")
    print("The shape of the water surface is the shape of the icosahedron's cross-section at this central, horizontal plane.")
    print("This plane cuts through the 6 faces that form the 'equator' of the icosahedron when it is oriented on a face.")
    print("The intersection of the plane with these 6 faces forms a regular six-sided polygon.")
    print("-" * 70)

    # Final Conclusion
    print("Conclusion: The shape of the water surface is a hexagon.")

# Execute the function to print the solution
solve_icosahedron_problem()