def explain_water_surface_shape():
    """
    Explains and determines the shape of the water surface in a half-filled icosahedron tank.
    """
    
    # Description of the object and setup
    object_name = "Icosahedron"
    base_shape = "equilateral triangle"
    fill_level = "half-filled"
    
    # Geometric Reasoning
    reasoning_step_1 = "An icosahedron is a centrally symmetric polyhedron. This means it has a center point, and any plane passing through this center divides its volume into two equal halves."
    reasoning_step_2 = "Since the tank is half-filled, the water surface must be a plane passing through the center of the icosahedron."
    reasoning_step_3 = "The tank is standing on one of its triangular faces. The water surface, under gravity, will be a horizontal plane parallel to this base."
    reasoning_step_4 = "The shape of the horizontal cross-section that passes through the center of an icosahedron oriented on a face is a regular hexagon."
    
    final_shape = "a regular hexagon"

    # Print the explanation
    print(f"To solve this problem, we follow these logical steps:")
    print(f"1. Object Properties: The tank is an {object_name}, which is a highly symmetric shape. It is standing on one of its {base_shape} faces.")
    print(f"2. Fill Condition: The tank is {fill_level}, meaning the water occupies half the total volume.")
    print(f"3. Symmetry: {reasoning_step_1}")
    print(f"4. Conclusion: {reasoning_step_2} When combining this with the fact that the water surface must be horizontal (parallel to the base), we find that the surface is the shape of the central horizontal cross-section of the icosahedron.")
    print(f"5. Final Shape: {reasoning_step_4}")
    print("-" * 20)
    print(f"Therefore, the shape of the water surface will be {final_shape}.")

explain_water_surface_shape()