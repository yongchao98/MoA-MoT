import math

def solve_cube_recycling_problem():
    """
    Calculates the number of chamfered cubes needed to recycle enough material for one new cube.
    """
    # Step 1: Define parameters and calculate the volume of a full cube.
    side_length = 10  # in mm
    volume_cube = side_length ** 3

    # Step 2: Calculate the volume of material removed from a single chamfer.
    # The chamfer is at 45 degrees, so its cross-section is a right-angled isosceles triangle.
    # The given sloping edge width of sqrt(2) is the hypotenuse of this triangle.
    # Let 'a' be the length of the triangle's other two equal sides.
    # From the Pythagorean theorem: a^2 + a^2 = (sqrt(2))^2  =>  2a^2 = 2  =>  a = 1 mm.
    chamfer_leg_length = 1.0  # in mm

    # The volume of the removed material (a triangular prism) is the cross-sectional area times the length.
    cross_section_area = 0.5 * chamfer_leg_length * chamfer_leg_length
    volume_per_chamfer = cross_section_area * side_length

    # Step 3: Calculate the total recycled volume from one cube.
    # Four edges on one face are chamfered.
    num_chamfered_edges = 4
    total_recycled_volume = num_chamfered_edges * volume_per_chamfer

    # Step 4: Determine the number of cubes needed.
    num_cubes_needed = volume_cube / total_recycled_volume

    # Final Output: Print the calculation and the result.
    print("The volume of one new cube is {} mm^3.".format(volume_cube))
    print("The volume of recycled material from one chamfered cube is {} mm^3.".format(int(total_recycled_volume)))
    print("The number of chamfered cubes needed to make one new cube is calculated as follows:")
    print("{} / {} = {}".format(int(volume_cube), int(total_recycled_volume), int(num_cubes_needed)))

solve_cube_recycling_problem()
<<<50>>>