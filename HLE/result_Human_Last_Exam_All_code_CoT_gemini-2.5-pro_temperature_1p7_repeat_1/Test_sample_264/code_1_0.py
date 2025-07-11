import math

def solve_intersection_volume():
    """
    Calculates the volume of the intersection of two cubes shifted along the main diagonal.
    """
    # Side length of the original cubes.
    L_cube = 1.0

    # The shift distance is 1/2.
    shift_distance = 0.5

    # The shift is along the main diagonal, whose direction vector is (1, 1, 1).
    # The length of this direction vector is sqrt(1^2 + 1^2 + 1^2) = sqrt(3).
    # To find the shift component 's' along each axis, we divide the shift distance
    # by the length of the diagonal vector.
    s = shift_distance / math.sqrt(3)

    # We can model the first cube as being centered at the origin, so it occupies the space
    # [-0.5, 0.5] on each axis. The second cube is shifted by (s, s, s).
    # The intersection of the two cubes is a smaller cube whose side length is given by:
    # L_intersect = L_cube - s
    L_intersect = L_cube - s

    # The volume of this intersection cube is its side length cubed.
    volume = L_intersect**3

    print("The volume of the intersection is calculated as (L_cube - s)^3.")
    print("Where L_cube = 1.0 (the side length of the cubes)")
    print(f"And s = {shift_distance} / sqrt(3) = {s}")
    print("\nThe final equation with the calculated values is:")
    print(f"({L_cube} - {s}) ** 3 = {volume}")

solve_intersection_volume()