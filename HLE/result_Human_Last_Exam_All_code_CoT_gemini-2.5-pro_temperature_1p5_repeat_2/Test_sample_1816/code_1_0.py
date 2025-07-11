import math

def solve_chamfered_cube_problem():
    """
    Calculates the number of chamfered cubes needed to recycle enough material
    to manufacture one new cube.
    """
    # Parameters from the problem description
    side_length = 10  # in mm
    num_chamfered_edges = 4
    chamfer_hypotenuse = math.sqrt(2)  # in mm

    # Step 1: Calculate the volume of one original cube.
    volume_cube = side_length ** 3

    # Step 2: Calculate the volume of the material removed by one chamfer.
    # The cross-section of the chamfer is an isosceles right-angled triangle.
    # The given width is the hypotenuse. Let the two equal sides be 'a'.
    # From Pythagoras theorem: a^2 + a^2 = hypotenuse^2
    # 2 * a^2 = (sqrt(2))^2 => 2 * a^2 = 2 => a = 1 mm.
    chamfer_side_a = math.sqrt(chamfer_hypotenuse**2 / 2)
    
    # The area of the triangular cross-section
    area_triangle = 0.5 * chamfer_side_a * chamfer_side_a
    
    # The volume of one chamfer cut (a triangular prism)
    volume_one_chamfer = area_triangle * side_length
    
    # The total recycled volume from one cube (since 4 edges are chamfered)
    volume_recycled_per_cube = num_chamfered_edges * volume_one_chamfer

    # Step 3: Calculate how many chamfered cubes are needed.
    # This is the volume of a new cube divided by the recycled volume per cube.
    # We use math.ceil to round up to the nearest whole number.
    num_cubes_needed = math.ceil(volume_cube / volume_recycled_per_cube)

    # Output the explanation and final answer
    print("This script calculates the number of chamfered cubes needed to make one new cube from recycled material.\n")
    print(f"The volume of one original cube is {side_length}mm * {side_length}mm * {side_length}mm = {volume_cube:.0f} mm^3.")
    print(f"The volume of material removed from one cube (by chamfering {num_chamfered_edges} edges) is {volume_recycled_per_cube:.0f} mm^3.")
    print("\nTo find the number of cubes needed, we perform the following calculation:")
    
    # The final equation with each number printed explicitly
    print(f"Number of cubes = (Volume of a new cube) / (Recycled volume per cube)")
    print(f"Number of cubes = {volume_cube:.0f} / {volume_recycled_per_cube:.0f} = {num_cubes_needed:.0f}")

solve_chamfered_cube_problem()
<<<50>>>