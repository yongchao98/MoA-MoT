import math

def solve_cube_recycling():
    """
    Calculates the number of chamfered cubes needed to recycle enough material
    to manufacture one new cube.
    """
    # Step 1: Define the initial parameters
    side_length = 10  # in mm
    num_chamfered_edges = 4
    chamfer_hypotenuse = math.sqrt(2)  # in mm

    # Step 2: Calculate the volume of the original cube
    volume_of_cube = side_length ** 3

    # Step 3: Calculate the volume of material removed by one chamfer
    # The cross-section of the chamfer is an isosceles right triangle.
    # Let 'a' be the length of the two equal sides.
    # a^2 + a^2 = chamfer_hypotenuse^2
    # 2 * a^2 = (sqrt(2))^2 = 2
    # a^2 = 1 => a = 1 mm
    chamfer_cut_side = 1.0  # in mm

    # The area of the triangular cross-section
    triangle_area = 0.5 * chamfer_cut_side * chamfer_cut_side

    # The volume of one chamfer (a triangular prism)
    volume_of_one_chamfer = triangle_area * side_length

    # Step 4: Calculate the total volume of recycled material from one cube
    total_recycled_volume_per_cube = volume_of_one_chamfer * num_chamfered_edges

    # Step 5: Calculate how many cubes are needed
    # This is the volume of a new cube divided by the recycled volume per cube
    num_cubes_needed = volume_of_cube / total_recycled_volume_per_cube

    print("The problem is to find how many chamfered cubes are needed to make a new one from the recycled material.")
    print(f"The volume of a new cube with side length {side_length} mm is {int(volume_of_cube)} mm^3.")
    print(f"The volume of material removed from one edge (one chamfer) is {int(volume_of_one_chamfer)} mm^3.")
    print(f"Since {num_chamfered_edges} edges are chamfered, the total recycled volume per cube is {int(total_recycled_volume_per_cube)} mm^3.")
    print("\nThe number of cubes needed is the volume of a new cube divided by the recycled volume per cube.")
    print("The final calculation is:")
    # The final instruction is to print the equation with the numbers
    print(f"{int(num_cubes_needed)} = {int(volume_of_cube)} / ({int(volume_of_one_chamfer)} * {num_chamfered_edges})")

solve_cube_recycling()
<<<50>>>