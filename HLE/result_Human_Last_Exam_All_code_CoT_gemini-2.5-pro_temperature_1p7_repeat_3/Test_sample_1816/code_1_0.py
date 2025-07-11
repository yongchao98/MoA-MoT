import math

def solve_chamfer_problem():
    """
    Calculates the number of chamfered cubes required to recycle
    enough material to manufacture one new cube.
    """

    # --- 1. Define constants and calculate the volume of one cube ---
    side_length = 10  # in mm
    volume_cube = side_length ** 3

    # --- 2. Calculate the volume of material removed from one edge ---
    # The chamfer is applied to 4 edges on one face.
    chamfered_edges_per_cube = 4
    # The cross-section of the cut is an isosceles right triangle.
    # The sloping edge (hypotenuse) width is sqrt(2).
    # Let 'a' be the length of the other two equal sides.
    # Using Pythagorean theorem: a^2 + a^2 = (sqrt(2))^2 => 2*a^2 = 2 => a = 1 mm.
    cut_depth_a = 1.0

    # The area of the triangular cross-section is (1/2) * base * height.
    area_cross_section = 0.5 * cut_depth_a * cut_depth_a

    # The length of the removed prism is the side length of the cube.
    volume_removed_per_edge = area_cross_section * side_length

    # --- 3. Calculate the total volume removed from one cube ---
    total_volume_removed_per_cube = volume_removed_per_edge * chamfered_edges_per_cube

    # --- 4. Calculate the number of cubes needed ---
    # N * total_volume_removed_per_cube = volume_cube
    num_cubes = volume_cube / total_volume_removed_per_cube

    # --- 5. Print the results and the final equation ---
    print(f"The volume of one new cube is {side_length}mm * {side_length}mm * {side_length}mm = {volume_cube} mm^3.")
    print(f"The volume of recycled material from one chamfered cube is {total_volume_removed_per_cube} mm^3.")
    print("To find the number of cubes needed, we set up the equation:")
    print(f"Number of Cubes * (Recycled Volume per Cube) = Volume of one New Cube")
    print("\nThe final equation is:")
    print(f"{int(num_cubes)} cubes = {int(volume_cube)} mm^3 / {int(total_volume_removed_per_cube)} mm^3 per cube")


solve_chamfer_problem()
<<<50>>>