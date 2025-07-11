import math

def solve_cube_recycling():
    """
    Calculates the number of chamfered cubes needed to recycle enough material
    to manufacture one new cube.
    """
    # Step 1: Define initial parameters
    side_length = 10  # in mm
    num_chamfered_edges = 4
    chamfer_hypotenuse = math.sqrt(2)  # in mm

    # Step 2: Calculate the volume of one original cube
    volume_of_cube = side_length ** 3
    print(f"The volume of one original cube is {side_length}^3 = {volume_of_cube} mm^3.")

    # Step 3: Calculate the volume of material removed from one chamfer
    # The cross-section of the cut is an isosceles right-angled triangle.
    # Let 'a' be the length of the two equal sides.
    # a^2 + a^2 = (sqrt(2))^2 => 2*a^2 = 2 => a^2 = 1 => a = 1 mm
    cut_side_length = math.sqrt(chamfer_hypotenuse**2 / 2)
    
    # The volume of one chamfer cut (a triangular prism)
    # Volume = Area of triangular base * length of the prism (cube's side)
    volume_of_one_chamfer = (0.5 * cut_side_length * cut_side_length) * side_length
    
    # Step 4: Calculate total recycled material from one cube
    total_recycled_volume = num_chamfered_edges * volume_of_one_chamfer
    print(f"The total volume of recycled material from one cube is {num_chamfered_edges} * {volume_of_one_chamfer} = {total_recycled_volume} mm^3.")

    # Step 5: Calculate the number of cubes needed
    # We need enough material for a new cube, so we divide the total volume needed
    # by the volume recycled per cube. We use math.ceil to ensure we have enough.
    num_cubes_needed = math.ceil(volume_of_cube / total_recycled_volume)

    print("\nTo find the number of cubes needed, we use the following equation:")
    print("Number of Cubes = Volume of one cube / Volume of recycled material from one cube")
    print(f"Therefore, the final calculation is:")
    print(f"{int(num_cubes_needed)} = {int(volume_of_cube)} / {int(total_recycled_volume)}")
    
    return num_cubes_needed

if __name__ == '__main__':
    answer = solve_cube_recycling()
    # The final answer is wrapped in <<<>>>
    print(f"\nThus, {int(answer)} chamfered cubes are needed to make one new cube.")
    print(f"<<<{int(answer)}>>>")
