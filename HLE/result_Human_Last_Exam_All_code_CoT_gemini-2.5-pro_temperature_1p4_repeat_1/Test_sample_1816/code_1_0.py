import math

def solve_chamfer_problem():
    """
    Calculates the number of chamfered cubes needed to recycle enough material for one new cube.
    """
    # Step 1: Define initial parameters and calculate the original cube's volume.
    cube_side_length = 10  # in mm
    volume_original_cube = cube_side_length ** 3

    # Step 2: Calculate the volume of a single chamfer cut.
    # The chamfer creates a sloping edge (hypotenuse of a triangle) with a width of sqrt(2).
    # The cross-section of the cut is an isosceles right-angled triangle, as the chamfer is at 45 degrees.
    # Let 'a' be the length of the two equal sides of the triangle.
    # By Pythagorean theorem: a^2 + a^2 = (sqrt(2))^2 => 2a^2 = 2 => a = 1 mm.
    chamfer_triangle_side = 1.0  # in mm

    # The area of the triangular cross-section of the cut.
    triangle_cross_section_area = 0.5 * chamfer_triangle_side * chamfer_triangle_side

    # The volume of the removed material for one chamfer is a triangular prism.
    # Its length is the side length of the cube.
    volume_one_chamfer = triangle_cross_section_area * cube_side_length

    # Step 3: Calculate the total recycled volume from one cube.
    # Four edges on one face are chamfered.
    num_chamfered_edges = 4
    total_removed_volume_per_cube = num_chamfered_edges * volume_one_chamfer

    # Step 4: Determine the number of cubes needed.
    # This is the original cube's volume divided by the recycled volume per cube.
    num_cubes_needed = volume_original_cube / total_removed_volume_per_cube

    # Print the final equation and the result.
    print("The volume of a new cube is {} mm^3.".format(int(volume_original_cube)))
    print("The recycled volume from one chamfered cube is {} mm^3.".format(int(total_removed_volume_per_cube)))
    print("To find the number of cubes needed, we perform the following calculation:")
    print("{} / {} = {}".format(int(volume_original_cube), int(total_removed_volume_per_cube), int(num_cubes_needed)))

solve_chamfer_problem()
<<<50>>>