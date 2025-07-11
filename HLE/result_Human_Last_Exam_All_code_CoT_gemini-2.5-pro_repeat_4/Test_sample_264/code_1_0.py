from fractions import Fraction

def solve_cube_intersection():
    """
    Calculates the volume of the intersection of two shifted cubes.
    """
    # Define the side length of the cubes as a fraction for precise output.
    side_length = Fraction(1, 1)

    # Define the shift along each axis. A shift "by 1/2 along the main diagonal"
    # is interpreted as a translation by the vector (1/2, 1/2, 1/2).
    shift = Fraction(1, 2)

    # Step 1: Define the coordinates of the two cubes.
    # Let the first cube (C1) occupy the region [0, 1] on each axis.
    c1_min = Fraction(0)
    c1_max = side_length

    # The second cube (C2) is shifted by `shift` from the first cube.
    c2_min = c1_min + shift
    c2_max = c1_max + shift

    print("Step 1: Define the cubes in a coordinate system.")
    print(f"Let the first cube be defined by the ranges: {c1_min} <= x,y,z <= {c1_max}.")
    print(f"The second cube is shifted by {shift} along each axis, so it is defined by the ranges: {c2_min} <= x,y,z <= {c2_max}.\n")

    # Step 2: Determine the dimensions of the intersection.
    # The intersection's range on each axis is the overlap of the two cubes' ranges.
    intersection_min = max(c1_min, c2_min)
    intersection_max = min(c1_max, c2_max)

    # The side length of the resulting intersection cube.
    intersection_side_length = intersection_max - intersection_min
    
    print("Step 2: Find the region of intersection.")
    print(f"The overlapping range for each axis is from max({c1_min}, {c2_min}) to min({c1_max}, {c2_max}), which is [{intersection_min}, {intersection_max}].")
    print(f"This means the intersection is a cube with a side length of {intersection_max} - {intersection_min} = {intersection_side_length}.\n")

    # Step 3: Calculate the volume of the intersection.
    volume = intersection_side_length ** 3

    print("Step 3: Calculate the volume of the intersection cube.")
    print("The volume is the side length cubed.")
    # The final equation with each number explicitly shown
    print(f"Volume = ({intersection_side_length}) * ({intersection_side_length}) * ({intersection_side_length}) = {volume}")

solve_cube_intersection()
<<<1/8>>>