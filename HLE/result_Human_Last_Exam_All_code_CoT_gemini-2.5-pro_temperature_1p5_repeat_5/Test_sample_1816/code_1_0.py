def solve_cube_recycling_problem():
    """
    Calculates the number of chamfered cubes needed to recycle enough
    material to manufacture one new cube.
    """
    # Define the initial parameters
    side_length = 10  # in mm
    num_chamfered_edges = 4

    # Step 1: Calculate the volume of a single cube
    volume_of_cube = side_length ** 3

    # Step 2: Calculate the volume of one chamfer
    # The chamfer is at 45 degrees, and the hypotenuse of the cut's
    # triangular cross-section is sqrt(2).
    # Using the Pythagorean theorem (a^2 + a^2 = (sqrt(2))^2), the length of the
    # other two sides of the right-angled isosceles triangle is 1 mm.
    cut_depth = 1  # in mm
    cross_section_area = 0.5 * cut_depth * cut_depth
    
    # The removed material is a triangular prism with length equal to the cube's side.
    volume_of_one_chamfer = cross_section_area * side_length

    # Step 3: Calculate total recycled material from one cube
    total_recycled_volume = num_chamfered_edges * volume_of_one_chamfer

    # Step 4: Calculate how many cubes are needed
    num_cubes_needed = volume_of_cube / total_recycled_volume
    
    # Print the explanation and the final equation with the calculated values
    print("The number of cubes needed is found by dividing the volume of a full cube by the total recycled volume from one chamfered cube.")
    print("\nCalculation steps:")
    print(f"1. The volume of a full cube with a side of {side_length} mm is {int(volume_of_cube)} mm^3.")
    print(f"2. The volume of material removed from one chamfer is {volume_of_one_chamfer} mm^3.")
    print(f"3. Since {num_chamfered_edges} edges are chamfered, the total recycled volume per cube is {num_chamfered_edges} * {volume_of_one_chamfer} = {total_recycled_volume} mm^3.")
    
    print("\nFinal Equation:")
    print(f"{int(volume_of_cube)} / ({num_chamfered_edges} * {int(volume_of_one_chamfer)}) = {int(num_cubes_needed)}")

solve_cube_recycling_problem()