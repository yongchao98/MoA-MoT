import math

def solve_cube_recycling():
    """
    Calculates the number of chamfered cubes needed to recycle enough material for one new cube.
    """
    # Parameters from the problem description
    side_length = 10  # mm
    num_chamfered_edges = 4
    chamfer_hypotenuse = math.sqrt(2)  # mm

    # Step 1: Calculate the volume of one original cube
    volume_cube = side_length ** 3

    # Step 2: Calculate the volume of the material removed from one edge (a triangular prism)
    # The base of the prism is a right-angled isosceles triangle.
    # The hypotenuse is given as sqrt(2). Let the equal sides be 'a'.
    # From Pythagorean theorem: a^2 + a^2 = (sqrt(2))^2 => 2*a^2 = 2 => a = 1 mm.
    cut_side_length = math.sqrt(chamfer_hypotenuse**2 / 2)

    # Area of the triangular base of the cut-off prism
    area_triangle = 0.5 * cut_side_length * cut_side_length

    # Volume of one cut-off prism (area of base * length of the prism)
    volume_per_chamfer = area_triangle * side_length

    # Step 3: Calculate the total recycled volume from one cube
    total_recycled_volume = num_chamfered_edges * volume_per_chamfer

    # Step 4: Calculate how many cubes are needed to make one new cube
    # This must be an integer, and since we need "enough" material, we
    # take the ceiling of the division.
    num_cubes_needed = math.ceil(volume_cube / total_recycled_volume)

    # Print the explanation and the final equation
    print("--- Calculation Steps ---")
    print(f"1. Volume of one original cube: {side_length}^3 = {volume_cube} mm^3")
    print(f"2. Volume of material removed per chamfered edge: (0.5 * {cut_side_length:.0f} * {cut_side_length:.0f}) * {side_length} = {volume_per_chamfer} mm^3")
    print(f"3. Total recycled material per cube: {num_chamfered_edges} * {volume_per_chamfer} = {total_recycled_volume} mm^3")
    print("\n--- Final Calculation ---")
    print("To find the number of cubes needed, we divide the volume of a new cube by the recycled volume per cube:")
    print(f"Number of Cubes = (Volume of one cube) / (Recycled material per cube)")
    # The final equation with all numbers
    print(f"Number of Cubes = {int(volume_cube)} / {int(total_recycled_volume)} = {int(num_cubes_needed)}")


solve_cube_recycling()
<<<50>>>