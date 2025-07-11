import math

def solve_cube_intersection_volume():
    """
    Calculates the volume of the intersection of two unit cubes
    shifted by 1/2 along the main diagonal.
    """
    # Step 1: Define the initial parameters.
    side_length = 1.0
    shift_distance = 0.5
    print("This script calculates the volume of the intersection of two shifted cubes.")
    print(f"The side length of each cube is L = {side_length}.")
    print(f"The shift distance along the main diagonal is d = {shift_distance}.")
    print("-" * 30)

    # Step 2: Calculate the properties of the main diagonal.
    # The main diagonal vector of a unit cube is (1, 1, 1).
    # Its length is sqrt(1^2 + 1^2 + 1^2) = sqrt(3).
    diagonal_length = math.sqrt(3)
    print("The shift is along the main diagonal, which has a direction vector of (1, 1, 1).")
    print(f"The length of this diagonal is sqrt(3) ≈ {diagonal_length:.6f}.")
    print("-" * 30)

    # Step 3: Calculate the shift component along each axis.
    # The total shift distance 'd' is distributed equally among the x, y, and z axes.
    # The shift component 's' along each axis is the total distance divided by the diagonal length.
    s = shift_distance / diagonal_length
    print("The shift component 's' along each axis is calculated as s = d / sqrt(3).")
    print(f"s = {shift_distance} / sqrt(3) ≈ {s:.6f}")
    print("-" * 30)

    # Step 4: Determine the side length of the intersection volume.
    # If the first cube is at [0, 1] on each axis, the shifted cube is at [s, 1+s].
    # The intersection on each axis is [s, 1].
    # The side length of the resulting intersection cube is (1 - s).
    intersection_side_length = side_length - s
    print("The intersection of the two cubes is a smaller cube.")
    print("Its side length is L_int = L - s.")
    print(f"L_int = {side_length} - {s:.6f} ≈ {intersection_side_length:.6f}")
    print("-" * 30)

    # Step 5: Calculate the final volume and print the equation.
    volume = intersection_side_length ** 3
    print("The volume of the intersection is V = (L - s)³.")
    print("\nFinal Equation:")
    print(f"Volume = ( {side_length} - ( {shift_distance} / sqrt(3) ) )³")
    print(f"Volume = ( {side_length} - {s:.6f} )³")
    print(f"Volume = ( {intersection_side_length:.6f} )³")
    print(f"Volume ≈ {volume:.6f}")

solve_cube_intersection_volume()
<<<0.359914>>>