import math

def solve_cube_recycling():
    """
    Calculates the number of chamfered cubes needed to recycle enough material
    for one new cube.
    """
    # 1. Define the initial parameters
    side_length = 10.0  # in mm
    num_chamfered_edges = 4

    # 2. Calculate the volume of a single uncut cube
    volume_cube = side_length ** 3
    
    # 3. Calculate the volume of material removed by one chamfer.
    # The chamfer is a 45-degree cut, creating a sloping edge (hypotenuse)
    # of sqrt(2) mm. This forms a right-angled isosceles triangle
    # as the cross-section of the cut.
    # Let 'x' be the length of the two equal sides of the triangle.
    # From Pythagorean theorem: x^2 + x^2 = (sqrt(2))^2 => 2*x^2 = 2 => x = 1 mm.
    chamfer_cut_side = 1.0
    
    # The volume of one chamfer cut is a triangular prism.
    # Volume = (Area of triangle base) * (length of the edge)
    area_triangle_base = 0.5 * chamfer_cut_side * chamfer_cut_side
    volume_one_chamfer = area_triangle_base * side_length

    # 4. Calculate the total volume of recycled material from one cube
    total_recycled_volume_per_cube = num_chamfered_edges * volume_one_chamfer
    
    # 5. Calculate how many chamfered cubes are needed
    # We need to round up because we can't have a fraction of a cube.
    num_cubes_needed = math.ceil(volume_cube / total_recycled_volume_per_cube)

    # Print the explanation and the final equation
    print(f"The volume of one full cube is {volume_cube:.0f} mm³.")
    print(f"The volume of material recycled from one chamfered cube is {total_recycled_volume_per_cube:.0f} mm³.")
    print(f"\nTo find the number of cubes needed, we divide the volume of a full cube by the recycled volume per cube:")
    print(f"Number of cubes = {volume_cube:.0f} / {total_recycled_volume_per_cube:.0f} = {num_cubes_needed:.0f}")
    
solve_cube_recycling()
<<<50>>>