import math

def solve_target_problem():
    """
    Calculates the expected number of hits on the yellow circles on a target.
    """
    # Define the given and derived parameters
    num_shots = 10000
    num_yellow_circles = 6
    r_yellow = 1
    
    # Based on the puzzle's geometric constraints, the dimensions are determined to be:
    target_width = 48
    target_height = 36

    # Calculate the areas
    total_yellow_area = num_yellow_circles * math.pi * (r_yellow ** 2)
    total_target_area = target_width * target_height

    # Calculate the expected number of hits
    expected_hits = num_shots * (total_yellow_area / total_target_area)

    # Print the explanation and the final equation with all its numbers
    print("To find the expected number of hits, we calculate the ratio of the yellow area to the total area and multiply by the number of shots.")
    print(f"\nNumber of shots: {num_shots}")
    print(f"Number of yellow circles: {num_yellow_circles}")
    print(f"Radius of one yellow circle: {r_yellow} cm")
    print(f"Total area of yellow circles = {num_yellow_circles} * pi * {r_yellow}^2 = {total_yellow_area:.4f} cm^2")
    print("\nTarget dimensions:")
    print(f"Width: {target_width} cm")
    print(f"Height: {target_height} cm")
    print(f"Total target area = {target_width} * {target_height} = {total_target_area} cm^2")
    
    print("\nThe final calculation is:")
    # Using formatted strings to display the equation with the numbers plugged in
    final_equation_str = f"{num_shots} * ({num_yellow_circles} * pi * {r_yellow}**2) / ({target_width} * {target_height})"
    print(final_equation_str)
    
    print(f"\nExpected hits = {expected_hits:.4f}")
    print(f"Rounding to the nearest whole number, we expect approximately {round(expected_hits)} hits.")

solve_target_problem()