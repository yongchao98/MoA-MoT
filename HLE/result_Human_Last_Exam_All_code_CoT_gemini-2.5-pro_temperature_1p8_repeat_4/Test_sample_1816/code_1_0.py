import math

def solve_cube_recycling_problem():
    """
    Calculates the number of chamfered cubes needed to recycle enough material for one new cube.
    """
    # Step 1: Define the properties of the cube and the chamfer process.
    side_length = 10  # mm
    chamfer_width = math.sqrt(2)  # mm, this is the hypotenuse of the cut triangle
    num_edges_chamfered = 4

    # Step 2: Calculate the volume of one original cube.
    # The volume of a cube is side_length^3.
    volume_of_one_cube = side_length ** 3
    print(f"The side length of the cube is {side_length} mm.")
    print(f"The volume of one original cube is {side_length} * {side_length} * {side_length} = {volume_of_one_cube} mm^3.")
    print("-" * 30)

    # Step 3: Calculate the volume of the material removed from one cube.
    # The chamfer cut has a cross-section of a right-angled isosceles triangle,
    # because the chamfer is at 45 degrees to the adjacent faces.
    # The chamfer width is the hypotenuse. Let 'a' be the length of the other two sides.
    # By Pythagorean theorem: a^2 + a^2 = hypotenuse^2
    # 2 * a^2 = (sqrt(2))^2 => 2 * a^2 = 2 => a = 1 mm.
    cut_depth = 1.0  # mm

    # The area of the triangular cross-section is (1/2) * base * height
    area_of_cut_triangle = 0.5 * cut_depth * cut_depth

    # The volume of material removed per edge is a triangular prism (Area * length)
    volume_removed_per_edge = area_of_cut_triangle * side_length

    # Total volume removed is for all four chamfered edges.
    total_volume_removed_per_cube = num_edges_chamfered * volume_removed_per_edge

    print(f"A chamfer is applied to {num_edges_chamfered} edges on one face.")
    print(f"The volume of material removed from a single edge (a triangular prism) is:")
    print(f"(0.5 * {cut_depth} * {cut_depth}) * {side_length} = {volume_removed_per_edge} mm^3.")
    print(f"The total volume of recycled material from one cube is:")
    print(f"{num_edges_chamfered} * {volume_removed_per_edge} = {total_volume_removed_per_cube} mm^3.")
    print("-" * 30)

    # Step 4: Calculate how many chamfered cubes are needed.
    # This is the total volume of one cube divided by the recycled volume per cube.
    num_cubes_needed = volume_of_one_cube / total_volume_removed_per_cube

    print("To find the number of chamfered cubes needed, we divide the volume of a full cube")
    print("by the recycled volume from one chamfered cube.")
    print("\nThe final equation is:")
    # Output the final equation with all the numbers, as requested.
    print(f"{int(num_cubes_needed)} = {int(volume_of_one_cube)} / {int(total_volume_removed_per_cube)}")

solve_cube_recycling_problem()
<<<50>>>