import math

def solve_cube_intersection():
    """
    Calculates the volume of the intersection of two cubes shifted along a main diagonal.
    """
    # Each cube has a side length of 1.
    side_length = 1

    # The shift along the main diagonal is 1/2.
    # In an axis-aligned system, this corresponds to a shift of 1/2 in each dimension.
    shift = 0.5

    # The intersection of two axis-aligned boxes is another axis-aligned box.
    # To find its volume, we first find the length of its sides.

    # Let the first cube occupy the range [0, 1] on an axis.
    # The second, shifted cube occupies the range [0.5, 1.5] on that same axis.
    # The intersection is the range [0.5, 1].
    # The length of this intersection is 1 - 0.5 = 0.5.
    
    # In general, for a side length S and a shift d, the overlap is S - d.
    intersection_side_length = side_length - shift

    # The volume of the intersection (which is a cube) is its side length cubed.
    volume = intersection_side_length ** 3

    print("Step 1: Define the initial parameters.")
    print(f"Side length of each cube = {side_length}")
    print(f"Shift along each axis = {shift}")
    print("\nStep 2: Determine the side length of the intersection volume.")
    print("The intersection forms a smaller cube. Its side length is the overlap in any single dimension.")
    print(f"Intersection side length = (Original side length) - (Shift)")
    print(f"Intersection side length = {side_length} - {shift} = {intersection_side_length}")
    print("\nStep 3: Calculate the volume of the intersection cube.")
    print("Volume = (Intersection side length)^3")
    print(f"Volume = ({side_length} - {shift})^3 = {intersection_side_length}^3 = {volume}")

solve_cube_intersection()