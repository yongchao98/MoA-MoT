def solve_icosahedron_problem():
    """
    This function explains the shape of the water surface in a half-filled icosahedron tank
    by reasoning through its geometry and symmetry.
    """
    print("Determining the shape of the water surface in a half-filled icosahedron tank:")
    print("=" * 75)

    # Step 1: Define the water level based on volume and symmetry.
    print("\nStep 1: The Water Level")
    print("An icosahedron resting on a face has a parallel top face. Due to its symmetry,")
    print("the horizontal plane that passes through its exact center divides the volume into two equal halves.")
    print("Therefore, when the tank is half-full, the water surface is this central, horizontal plane.")

    # Step 2: Describe the cross-section created by the water surface.
    print("\nStep 2: The Cross-Section")
    print("The question is now: What is the shape of the icosahedron's cross-section at its center,")
    print("parallel to a face?")
    
    # Step 3: Count the sides of the resulting shape.
    num_sides = 6
    print(f"\nStep 3: Finding the Number of Sides")
    print(f"The central plane intersects the edges that connect the top half of the icosahedron to the bottom half.")
    print(f"By examining the structure of the icosahedron, we find there are exactly {num_sides} such edges.")
    print(f"The intersection points on these {num_sides} edges form the vertices of the water surface.")
    
    # Step 4: Conclude the shape.
    print("\nStep 4: The Final Shape")
    print(f"A polygon with {num_sides} vertices and sides is a hexagon.")
    print("Because of the perfect 3-fold rotational and central symmetry of the icosahedron,")
    print("this polygon must be a regular hexagon.")
    print("\n---------------------------------------------------------------------------")
    print("Conclusion: The shape of the water surface is a regular hexagon.")
    print("---------------------------------------------------------------------------")


# Run the reasoning process.
solve_icosahedron_problem()