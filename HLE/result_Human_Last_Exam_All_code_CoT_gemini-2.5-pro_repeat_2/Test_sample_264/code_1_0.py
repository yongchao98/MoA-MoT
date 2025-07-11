import math

def solve_intersection_volume():
    """
    Calculates the volume of the intersection of two cubes shifted along the main diagonal.
    """
    # Side length of each cube
    cube_side_length = 1.0

    # The magnitude of the shift along the main diagonal
    shift_magnitude = 0.5

    # The main diagonal is in the direction (1, 1, 1). The length of this direction vector is sqrt(3).
    # We need to find the component 'd' of the shift vector along each axis.
    # The shift vector s = (d, d, d), and its magnitude is ||s|| = sqrt(d^2 + d^2 + d^2) = d*sqrt(3).
    # We are given ||s|| = 0.5.
    # So, d*sqrt(3) = 0.5, which means d = 0.5 / sqrt(3).
    d = shift_magnitude / math.sqrt(3)

    # The intersection of the two cubes is a smaller cube.
    # Its side length is the original side length minus the shift component 'd'.
    intersect_side_length = cube_side_length - d

    # The volume of the intersection is the cube of its side length.
    volume = intersect_side_length ** 3

    print("Step 1: Define the initial parameters.")
    print(f"Side length of the cubes (L): {cube_side_length}")
    print(f"Shift magnitude along the main diagonal: {shift_magnitude}\n")

    print("Step 2: Calculate the shift component 'd' along each axis.")
    print(f"d = shift_magnitude / sqrt(3) = {shift_magnitude} / {math.sqrt(3):.6f} = {d:.6f}\n")

    print("Step 3: Calculate the side length of the intersection volume.")
    print(f"Intersection side length = L - d = {cube_side_length} - {d:.6f} = {intersect_side_length:.6f}\n")

    print("Step 4: Calculate the final volume.")
    # Show the equation with numbers
    print("The final equation for the volume is (L - d)^3.")
    print(f"Volume = ({cube_side_length} - {shift_magnitude} / sqrt(3))^3")
    print(f"Volume = ({cube_side_length:.1f} - {d:.6f})^3 = ({intersect_side_length:.6f})^3")
    print(f"Final Volume: {volume:.6f}")

solve_intersection_volume()