import math

def solve_cube_recycling():
    """
    Calculates the number of chamfered cubes needed to recycle enough material
    to manufacture one new cube.
    """
    # 1. Define initial parameters
    side_length = 10  # in mm
    num_chamfered_edges = 4
    chamfer_hypotenuse = math.sqrt(2)  # in mm

    # 2. Calculate the volume of one original cube
    volume_cube = side_length ** 3
    print(f"The volume of one original cube is {side_length}mm * {side_length}mm * {side_length}mm = {volume_cube} mm^3.")

    # 3. Calculate the dimensions and volume of the material cut from one edge.
    # For a 45-degree right triangle with hypotenuse h, the two other sides s are equal.
    # s^2 + s^2 = h^2  => 2*s^2 = h^2  => s = sqrt(h^2 / 2)
    cut_side_length_sq = (chamfer_hypotenuse ** 2) / 2
    cut_side_length = math.sqrt(cut_side_length_sq)

    # The cross-section is a right triangle. Area = 0.5 * base * height
    cross_section_area = 0.5 * cut_side_length * cut_side_length

    # Volume of the cut-off prism = area * length of the edge
    volume_one_chamfer = cross_section_area * side_length
    print(f"The volume of material removed from one chamfered edge is {volume_one_chamfer} mm^3.")

    # 4. Calculate the total recycled volume from one cube
    total_recycled_volume = num_chamfered_edges * volume_one_chamfer
    print(f"The total recycled volume from one cube (with {num_chamfered_edges} edges chamfered) is {total_recycled_volume} mm^3.")

    # 5. Calculate how many cubes are needed
    num_cubes_needed = volume_cube / total_recycled_volume
    
    print("\nTo find the number of cubes needed, we use the following equation:")
    print("Number of Cubes = (Volume of one Cube) / (Volume from one Chamfer * Number of Chamfered Edges)")
    print(f"Number of Cubes = {int(volume_cube)} / ({int(volume_one_chamfer)} * {int(num_chamfered_edges)})")
    print(f"Number of Cubes = {int(volume_cube)} / {int(total_recycled_volume)}")
    print(f"Therefore, the final number of cubes needed is: {int(num_cubes_needed)}")

solve_cube_recycling()