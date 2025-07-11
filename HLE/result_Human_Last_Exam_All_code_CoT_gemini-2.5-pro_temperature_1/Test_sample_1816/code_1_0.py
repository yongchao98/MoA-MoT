import math

def solve_chamfered_cube_problem():
    """
    Calculates the number of chamfered cubes needed to recycle enough material
    for one new cube.
    """
    # Step 1: Define initial parameters and calculate the volume of a single cube.
    side_length = 10  # in mm
    chamfer_width = math.sqrt(2) # in mm
    num_chamfered_edges = 4

    volume_cube = side_length ** 3
    print(f"The volume of a single solid cube is {side_length}mm * {side_length}mm * {side_length}mm = {volume_cube} mm^3.")
    print("-" * 20)

    # Step 2: Analyze the chamfer geometry.
    # The cross-section is a right-angled isosceles triangle.
    # chamfer_width is the hypotenuse (w). The legs are 'a'.
    # a^2 + a^2 = w^2  =>  2a^2 = (sqrt(2))^2  =>  2a^2 = 2  => a = 1
    chamfer_leg_a = math.sqrt(chamfer_width**2 / 2)
    cross_section_area = 0.5 * chamfer_leg_a * chamfer_leg_a
    print(f"The chamfer has a cross-section of a right-angled isosceles triangle.")
    print(f"The leg of this triangle ('a') is {chamfer_leg_a:.1f} mm.")
    print(f"The cross-sectional area of the cut is 0.5 * {chamfer_leg_a:.1f} * {chamfer_leg_a:.1f} = {cross_section_area:.1f} mm^2.")
    print("-" * 20)
    
    # Step 3: Calculate the volume of removed material per cube.
    # We use the principle of inclusion-exclusion: V_total = V_sum - V_overlap
    
    # First, calculate the sum of volumes of 4 ideal prisms without considering overlap.
    volume_one_prism = cross_section_area * side_length
    total_prism_volume = num_chamfered_edges * volume_one_prism
    print(f"The volume of one ideal chamfer prism is {cross_section_area:.1f} mm^2 * {side_length} mm = {volume_one_prism} mm^3.")
    print(f"The sum of volumes for {num_chamfered_edges} ideal prisms is {num_chamfered_edges} * {volume_one_prism} = {total_prism_volume} mm^3.")
    
    # Second, calculate the volume of the overlap at one corner.
    # The overlap volume is the integral of (a-z)^2 dz from 0 to a, which is a^3 / 3.
    # Here a = chamfer_leg_a = 1.
    volume_one_overlap = (chamfer_leg_a**3) / 3
    total_overlap_volume = num_chamfered_edges * volume_one_overlap
    print(f"The volume of the overlap at each of the {num_chamfered_edges} corners is {volume_one_overlap:.4f} mm^3.")
    
    # Finally, calculate the actual volume removed per cube.
    volume_removed_per_cube = total_prism_volume - total_overlap_volume
    print(f"\nThe actual volume of recycled material from one cube is {total_prism_volume} - {num_chamfered_edges} * {volume_one_overlap:.4f} = {volume_removed_per_cube:.4f} mm^3.")
    print("-" * 20)

    # Step 4: Determine the number of cubes needed.
    # N * volume_removed_per_cube >= volume_cube
    num_cubes_float = volume_cube / volume_removed_per_cube
    num_cubes_required = math.ceil(num_cubes_float)
    
    print("To find the number of cubes (N) needed, we solve the equation:")
    print(f"N * {volume_removed_per_cube:.4f} mm^3 >= {volume_cube} mm^3")
    print(f"N >= {volume_cube} / {volume_removed_per_cube:.4f}")
    print(f"N >= {num_cubes_float:.4f}")
    print("\nSince the number of cubes must be an integer, we round up.")
    print(f"The minimum number of chamfered cubes required is {num_cubes_required}.")

    # Final equation check
    total_recycled_volume = num_cubes_required * volume_removed_per_cube
    print(f"\nFinal check: {num_cubes_required} cubes * {volume_removed_per_cube:.4f} mm^3/cube = {total_recycled_volume:.2f} mm^3, which is >= {volume_cube} mm^3.")
    
    return num_cubes_required

if __name__ == '__main__':
    answer = solve_chamfered_cube_problem()
    # The final answer is returned in the specified format.
    print(f"\n<<<{answer}>>>")