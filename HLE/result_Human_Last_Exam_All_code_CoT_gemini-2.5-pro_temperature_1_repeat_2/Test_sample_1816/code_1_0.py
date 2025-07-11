import math

def solve_cube_recycling_problem():
    """
    Calculates the number of chamfered cubes needed to recycle enough material
    to manufacture one new cube.
    """
    # Step 1: Define initial parameters and calculate the volume of the original cube.
    side_length = 10.0  # in mm
    volume_cube = side_length ** 3
    print(f"The volume of one original cube is {side_length}^3 = {volume_cube} mm^3.")

    # Step 2: Calculate the volume of material removed from one edge.
    # The chamfer has a width of sqrt(2), which is the hypotenuse of the
    # triangular cross-section. The cross-section is an isosceles right triangle.
    # Let the two equal sides be 'x'. From Pythagoras: x^2 + x^2 = (sqrt(2))^2 => 2*x^2 = 2 => x = 1.
    chamfer_hypotenuse = math.sqrt(2)
    cut_depth = math.sqrt(chamfer_hypotenuse**2 / 2.0)
    
    # The removed shape is a triangular prism.
    # Area of the triangular base = 0.5 * base * height
    area_triangle = 0.5 * cut_depth * cut_depth
    # Volume of the prism = area of base * length
    volume_one_chamfer = area_triangle * side_length
    
    # Step 3: Calculate the total volume of material removed from one cube.
    # The process is applied to 4 edges.
    num_chamfered_edges = 4
    total_volume_removed = volume_one_chamfer * num_chamfered_edges
    print(f"The volume of recycled material from one cube is {num_chamfered_edges} * {volume_one_chamfer} = {total_volume_removed} mm^3.")

    # Step 4: Calculate how many cubes are needed.
    # This is the volume of a full cube divided by the recycled volume per cube.
    num_cubes_needed = volume_cube / total_volume_removed

    print("\nTo find the number of cubes needed, we use the following equation:")
    print(f"Number of Cubes = (Volume of a new cube) / (Recycled volume per cube)")
    # Outputting each number in the final equation as requested
    print(f"Number of Cubes = {int(volume_cube)} / {int(total_volume_removed)}")
    
    print(f"\nTherefore, {int(num_cubes_needed)} chamfered cubes are needed.")

solve_cube_recycling_problem()

<<<50>>>