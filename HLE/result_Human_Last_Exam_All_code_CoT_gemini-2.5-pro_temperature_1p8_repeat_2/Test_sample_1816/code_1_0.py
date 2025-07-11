import math

def solve_cube_recycling_problem():
    """
    This function calculates the number of chamfered cubes needed to recycle
    enough material for one new cube.
    """

    # Step 1: Define the initial parameters of the cube in mm.
    side_length = 10.0
    chamfered_edges_per_cube = 4.0
    
    # Step 2: Calculate the volume of a single cube.
    # Formula: V = side_length^3
    volume_of_one_cube = side_length ** 3
    
    # Step 3: Calculate the volume of material removed from a single edge (a triangular prism).
    # The cross-section of the cut is a right-angled triangle.
    # The hypotenuse is sqrt(2) mm. Since the chamfer is at 45 degrees, the other two sides are equal.
    # Let the side be 'a'. a^2 + a^2 = (sqrt(2))^2 => 2*a^2 = 2 => a = 1 mm.
    cut_depth = 1.0
    
    # The area of the triangular cross-section is 0.5 * base * height
    area_of_triangular_cut = 0.5 * cut_depth * cut_depth
    
    # The volume of one cut is the area of the triangle * the length of the edge.
    volume_of_one_chamfer = area_of_triangular_cut * side_length
    
    # Step 4: Calculate the total recycled volume from one chamfered cube.
    total_recycled_volume_per_cube = volume_of_one_chamfer * chamfered_edges_per_cube
    
    # Step 5: Calculate how many chamfered cubes are needed to make one new cube.
    number_of_cubes_needed = volume_of_one_cube / total_recycled_volume_per_cube

    # Step 6: Print the final calculation and answer, showing each number.
    print(f"The volume of one new cube is {int(volume_of_one_cube)} mm^3.")
    print(f"The volume of recycled material from one chamfered cube is {int(total_recycled_volume_per_cube)} mm^3.")
    print("The number of chamfered cubes needed is calculated by dividing the volume of a new cube by the recycled volume per cube.")
    print("Final Calculation:")
    print(f"{int(volume_of_one_cube)} / ({int(chamfered_edges_per_cube)} * {int(volume_of_one_chamfer)}) = {int(number_of_cubes_needed)}")

solve_cube_recycling_problem()
<<<50>>>