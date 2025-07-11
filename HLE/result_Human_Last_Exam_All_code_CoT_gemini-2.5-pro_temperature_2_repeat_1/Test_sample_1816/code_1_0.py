import math

def solve_cube_recycling():
    """
    Calculates the number of chamfered cubes needed to recycle enough material
    to manufacture one new cube.
    """

    # Step 1: Define parameters and calculate the volume of the original cube.
    side_length = 10  # in mm
    volume_of_one_cube = side_length ** 3
    print(f"The volume of one original cube is {side_length}^3 = {volume_of_one_cube} mm^3.")

    # Step 2: Calculate the volume of material removed from a single chamfered edge.
    # The chamfer is at 45 degrees, forming a cut with a right-angled isosceles triangular cross-section.
    # The given sloping edge width (sqrt(2) mm) is the hypotenuse of this triangle.
    # Let 'a' be the length of the two equal sides of the triangle cut from the faces.
    # From Pythagoras' theorem: a^2 + a^2 = (sqrt(2))^2  =>  2a^2 = 2  =>  a = 1 mm.
    a = 1.0
    
    # The area of the triangular cross-section is (1/2) * base * height
    area_of_cut_triangle = 0.5 * a * a
    
    # The volume of the removed prism is the cross-section area times the length of the cube's edge.
    volume_of_one_cut = area_of_cut_triangle * side_length
    print(f"The volume of material removed from one edge is {volume_of_one_cut} mm^3.")

    # Step 3: Calculate the total recycled volume from one cube.
    # Four edges on one face are chamfered.
    num_chamfered_edges = 4
    total_recycled_volume = num_chamfered_edges * volume_of_one_cut
    print(f"The total recycled material from one cube (with {num_chamfered_edges} chamfers) is {num_chamfered_edges} * {volume_of_one_cut} = {total_recycled_volume} mm^3.")

    # Step 4: Calculate the number of cubes needed.
    # This is the volume of a full cube divided by the recycled volume per cube.
    num_cubes_needed = volume_of_one_cube / total_recycled_volume
    
    print("\nTo find the number of chamfered cubes needed to make a new cube, we perform the following calculation:")
    print(f"Number of Cubes = (Volume of one cube) / (Recycled material per cube)")
    print(f"Number of Cubes = {int(volume_of_one_cube)} / {int(total_recycled_volume)} = {int(num_cubes_needed)}")


solve_cube_recycling()
<<<50>>>