import math

def solve_cube_recycling():
    """
    Calculates the number of chamfered cubes needed to recycle enough material
    to manufacture one new cube.
    """
    # Step 1: Define the properties of the cube and the chamfer
    cube_side_length = 10  # in mm
    num_chamfered_edges = 4
    chamfer_hypotenuse = math.sqrt(2)  # in mm

    # Step 2: Calculate the volume of one original cube
    volume_full_cube = cube_side_length ** 3

    # Step 3: Calculate the volume of material removed from one edge (one chamfer)
    # The cross-section of the chamfer cut is an isosceles right-angled triangle.
    # Let 'a' be the length of the two equal sides.
    # a^2 + a^2 = hypotenuse^2  => 2 * a^2 = (sqrt(2))^2 => 2 * a^2 = 2 => a^2 = 1 => a = 1 mm
    chamfer_leg_a = math.sqrt(chamfer_hypotenuse**2 / 2)
    
    # Area of the triangular cross-section = 0.5 * base * height
    triangle_area = 0.5 * chamfer_leg_a * chamfer_leg_a
    
    # Volume of one chamfer (a triangular prism) = Area * length
    volume_one_chamfer = triangle_area * cube_side_length
    
    # Step 4: Calculate the total recycled volume from one chamfered cube
    total_recycled_volume_per_cube = volume_one_chamfer * num_chamfered_edges
    
    # Step 5: Calculate how many cubes are needed
    # This is the total volume needed divided by the recycled volume per cube.
    # We use math.ceil because we need enough material for a *whole* new cube.
    num_cubes_needed = math.ceil(volume_full_cube / total_recycled_volume_per_cube)

    # Output the explanation and the final equation
    print("--- Calculation ---")
    print(f"Volume of one original cube: {cube_side_length} * {cube_side_length} * {cube_side_length} = {volume_full_cube} mm^3")
    print(f"Volume recycled from one cube: {num_chamfered_edges} edges * ({triangle_area} mm^2 * {cube_side_length} mm/edge) = {total_recycled_volume_per_cube} mm^3")
    print("\n--- Final Equation ---")
    print("Number of cubes needed = (Volume of a full cube) / (Recycled volume per cube)")
    print(f"Number of cubes needed = {int(volume_full_cube)} / {int(total_recycled_volume_per_cube)} = {int(num_cubes_needed)}")
    print("\nTherefore, the number of chamfered cubes needed is:")
    print(int(num_cubes_needed))

solve_cube_recycling()
<<<50>>>