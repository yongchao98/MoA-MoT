import math

def solve_cube_recycling_problem():
    """
    Calculates the number of chamfered cubes needed to recycle enough material
    for one new cube.
    """

    # Problem parameters
    side_length = 10  # in mm
    num_chamfered_edges = 4
    # The width of the sloping edge is the hypotenuse of the triangular cross-section
    chamfer_hypotenuse = math.sqrt(2)

    # Step 1: Calculate the volume of a single, uncut cube
    volume_cube = side_length ** 3
    print("Step 1: Calculate the volume of one original cube.")
    print(f"The volume of a cube with a side length of {side_length} mm is:")
    print(f"V_cube = {side_length} * {side_length} * {side_length} = {volume_cube} mm^3\n")

    # Step 2: Calculate the volume of material removed from one cube
    # The chamfer has a 45-degree angle, so its cross-section is an isosceles
    # right-angled triangle. The chamfer_hypotenuse is the hypotenuse.
    # Let 'a' be the length of the two equal sides of the triangle.
    # From Pythagoras theorem: a^2 + a^2 = hypotenuse^2
    # 2*a^2 = (sqrt(2))^2 => 2*a^2 = 2 => a^2 = 1 => a = 1 mm
    cut_depth_a = math.sqrt(chamfer_hypotenuse**2 / 2)
    
    # The volume of one chamfer cut is a triangular prism
    # Volume = Area of triangular base * length
    volume_one_chamfer = (0.5 * cut_depth_a * cut_depth_a) * side_length
    
    # Total recycled volume from one cube is from all chamfered edges
    total_recycled_volume = num_chamfered_edges * volume_one_chamfer
    
    print("Step 2: Calculate the recycled volume from one chamfered cube.")
    print(f"The material removed by one chamfer is a triangular prism {side_length} mm long.")
    print(f"Its cross-section is a right triangle with a hypotenuse of {chamfer_hypotenuse:.4f} mm.")
    print(f"The other two sides of the triangle are each {cut_depth_a} mm long.")
    print(f"Volume from one chamfer = (1/2 * {cut_depth_a} * {cut_depth_a}) * {side_length} = {volume_one_chamfer} mm^3")
    print(f"Since {num_chamfered_edges} edges are chamfered, the total recycled volume per cube is:")
    print(f"V_recycled = {num_chamfered_edges} * {volume_one_chamfer} = {total_recycled_volume} mm^3\n")
    
    # Step 3: Calculate the number of cubes needed
    # num_cubes_needed = volume_cube / total_recycled_volume
    # We use math.ceil because we need a whole number of cubes.
    num_cubes_needed = math.ceil(volume_cube / total_recycled_volume)

    print("Step 3: Calculate the number of cubes needed.")
    print("This is the volume of a new cube divided by the recycled volume per cube.")
    print("Number of cubes = (Volume of new cube) / (Recycled volume per cube)")
    # The final equation with numbers
    print(f"Number of cubes = {int(volume_cube)} / {int(total_recycled_volume)}")
    print(f"Number of cubes = {int(num_cubes_needed)}\n")

    print(f"Therefore, {int(num_cubes_needed)} chamfered cubes are needed to make another cube.")

solve_cube_recycling_problem()
<<<50>>>