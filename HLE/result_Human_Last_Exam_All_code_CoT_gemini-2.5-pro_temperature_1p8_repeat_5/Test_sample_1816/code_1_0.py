import math

def solve_chamfer_problem():
    """
    Calculates the number of chamfered cubes needed to recycle enough material
    to manufacture one new cube.
    """
    # --- Parameters ---
    side_length = 10  # in mm
    num_chamfered_edges = 4
    # The width of the sloping edge (hypotenuse of the triangular cross-section)
    chamfer_hypotenuse = math.sqrt(2)

    # --- Step 1: Calculate the volume of a full cube ---
    volume_cube = side_length ** 3

    # --- Step 2: Calculate the volume of a single chamfered piece ---
    # The cross-section is a right-angled isosceles triangle.
    # From the Pythagorean theorem: a^2 + a^2 = hypotenuse^2
    # 2*a^2 = (sqrt(2))^2 => 2*a^2 = 2 => a^2 = 1 => a = 1 mm.
    # 'a' is the side length of the triangular cross-section cut from the face.
    cut_depth = math.sqrt(chamfer_hypotenuse**2 / 2)

    # The area of the triangular cross-section is 0.5 * base * height
    area_triangle = 0.5 * cut_depth * cut_depth

    # The volume of the removed triangular prism is its area times its length
    volume_one_chamfer = area_triangle * side_length

    # --- Step 3: Calculate the total recycled volume from one cube ---
    total_recycled_volume_per_cube = num_chamfered_edges * volume_one_chamfer

    # --- Step 4: Calculate the number of cubes needed ---
    # We divide the volume of a full cube by the recycled material per cube.
    # math.ceil is used to ensure we have *enough* material, rounding up if necessary.
    num_cubes_needed = math.ceil(volume_cube / total_recycled_volume_per_cube)

    # --- Print the results ---
    print(f"The volume of one original cube is {volume_cube} mm^3.")
    print(f"The volume of recycled material from one chamfered edge is {volume_one_chamfer} mm^3.")
    print(f"The total recycled volume from one cube with {num_chamfered_edges} chamfers is {total_recycled_volume_per_cube} mm^3.")
    print("\nTo find how many cubes are needed, we perform the following calculation:")
    # The final equation with each number printed explicitly
    print(f"Volume of a full cube / (Number of chamfers * Volume per chamfer)")
    print(f"{int(volume_cube)} / ({num_chamfered_edges} * {int(volume_one_chamfer)}) = {int(num_cubes_needed)}")
    print(f"\nTherefore, {int(num_cubes_needed)} chamfered cubes are needed to accumulate enough material for a new cube.")

if __name__ == "__main__":
    solve_chamfer_problem()

<<<50>>>