import math

def calculate_intersection_volume():
    """
    Calculates and explains the volume of the intersection of two shifted unit cubes.
    """
    # Define problem parameters
    cube_side_length = 1.0
    shift_magnitude = 0.5

    print("Step 1: Understanding the geometry")
    print(f"Let the first cube C1 be defined by 0 <= x, y, z <= {cube_side_length}.")
    print("The main diagonal connects (0,0,0) to (1,1,1), so its direction vector is (1,1,1).")
    print("-" * 30)

    # Calculate the shift component along each axis
    # The length of the vector (1,1,1) is sqrt(1^2 + 1^2 + 1^2) = sqrt(3)
    diagonal_vector_length = math.sqrt(3)
    shift_along_axis = shift_magnitude / diagonal_vector_length

    print("Step 2: Calculating the shift per axis")
    print(f"The shift has a magnitude of {shift_magnitude} along the main diagonal.")
    print("The shift component along each axis is the magnitude divided by the length of the diagonal vector (sqrt(3)).")
    print(f"Shift along one axis (d) = {shift_magnitude} / sqrt(3) = {shift_along_axis}")
    print("-" * 30)

    # Determine the intersection's side length
    # The first cube C1 is in the range [0, 1].
    # The second cube C2 is shifted by 'd', so its range is [d, 1+d].
    # The intersection range is [max(0,d), min(1, 1+d)] which is [d, 1].
    # The side length of the intersection is (1 - d).
    intersection_side_length = cube_side_length - shift_along_axis
    
    print("Step 3: Determining the intersection dimensions")
    print("The intersection of the two cubes is a smaller cube.")
    print(f"Its side length is the original side length minus the shift along one axis.")
    print(f"Intersection side length (L) = {cube_side_length} - {shift_along_axis} = {intersection_side_length}")
    print("-" * 30)
    
    # Calculate the final volume
    intersection_volume = intersection_side_length ** 3

    print("Step 4: Calculating the final volume")
    print("The volume of the intersection is its side length cubed (V = L^3).")
    print("\nFinal Equation:")
    print(f"Volume = ({cube_side_length} - {shift_along_axis})^3")
    print(f"Volume = ({intersection_side_length})^3")
    print(f"Volume = {intersection_volume}")

# Run the calculation
calculate_intersection_volume()