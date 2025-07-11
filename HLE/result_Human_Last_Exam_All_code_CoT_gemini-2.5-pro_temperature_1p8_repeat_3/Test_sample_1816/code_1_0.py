import math

def solve_cube_recycling():
    """
    Calculates the number of chamfered cubes needed to recycle enough material for one new cube.
    """
    # Step 1: Define the properties of the cube and the chamfer
    side_length = 10.0  # mm
    num_chamfered_edges = 4
    chamfer_hypotenuse = math.sqrt(2) # mm

    # Step 2: Calculate the volume of a single, original cube
    volume_cube = side_length ** 3

    # Step 3: Calculate the volume of material removed by one chamfer
    # The chamfer is at 45 degrees, so its cross-section is a right-angled isosceles triangle.
    # Let 'a' be the length of the two equal sides of the triangle.
    # By Pythagorean theorem: a^2 + a^2 = hypotenuse^2  =>  2*a^2 = (sqrt(2))^2 = 2
    # This means a^2 = 1, so the side 'a' is 1 mm.
    chamfer_leg_length = math.sqrt(chamfer_hypotenuse**2 / 2.0)
    
    # The area of the triangular cross-section of the cut
    area_triangle = 0.5 * chamfer_leg_length * chamfer_leg_length
    
    # The volume of one removed prism (chamfer) = Area * length of the edge
    volume_one_chamfer = area_triangle * side_length
    
    # Step 4: Calculate the total recycled volume from one cube
    total_volume_removed_per_cube = volume_one_chamfer * num_chamfered_edges
    
    # Step 5: Calculate how many chamfered cubes are needed to make one new cube
    num_cubes_needed = volume_cube / total_volume_removed_per_cube
    
    # Step 6: Print the final equation and result
    print("To find the number of cubes needed, we use the following equation:")
    print("Number of Cubes = (Volume of one Cube) / (Volume recycled from one Cube)")
    print("\nBreaking it down:")
    print(f"Volume of one Cube = {side_length:.0f} * {side_length:.0f} * {side_length:.0f} = {volume_cube:.0f} mm^3")
    print(f"Volume recycled from one Cube = {num_chamfered_edges} * (Area of chamfer's triangle * length)")
    print(f"Volume recycled from one Cube = {num_chamfered_edges} * ( (0.5 * {chamfer_leg_length:.0f} * {chamfer_leg_length:.0f}) * {side_length:.0f} ) = {total_volume_removed_per_cube:.0f} mm^3")
    print("\nFinal Calculation:")
    print(f"Number of Cubes = {volume_cube:.0f} / {total_volume_removed_per_cube:.0f} = {num_cubes_needed:.0f}")

solve_cube_recycling()
<<<50>>>