import math

def solve_prime_path():
    """
    Calculates the number of distinct 4-move paths from (1, 1) to (5, 7)
    on the PrimeGrid+1.
    """
    # The set of allowed coordinate values
    allowed_coords = {1, 2, 3, 5, 7}

    # --- Case 1: 1 Horizontal move and 3 Vertical moves ---
    
    # We need to count the number of 3-move vertical paths from y_start=1 to y_end=7.
    # The path sequence is y_start -> y1 -> y2 -> y_end.
    # Constraints:
    # 1. y1 must be different from y_start.
    # 2. y2 must be different from y1.
    # 3. y_end must be different from y2.
    y_start, y_end = 1, 7
    v_paths_count = 0
    # First intermediate y-coordinate
    for y1 in allowed_coords - {y_start}:
        # Second intermediate y-coordinate
        for y2 in allowed_coords - {y1, y_end}:
            v_paths_count += 1

    # The single horizontal move can be placed in 4 different positions in the sequence.
    num_placements_for_h = 4
    total_case1 = num_placements_for_h * v_paths_count

    # --- Case 2: 3 Horizontal moves and 1 Vertical move ---

    # We need to count the number of 3-move horizontal paths from x_start=1 to x_end=5.
    # The path sequence is x_start -> x1 -> x2 -> x_end.
    x_start, x_end = 5, 7  # The destination is (5, 7)
    x_start = 1
    x_end = 5
    h_paths_count = 0
    # First intermediate x-coordinate
    for x1 in allowed_coords - {x_start}:
        # Second intermediate x-coordinate
        for x2 in allowed_coords - {x1, x_end}:
            h_paths_count += 1
            
    # The single vertical move can be placed in 4 different positions.
    num_placements_for_v = 4
    total_case2 = num_placements_for_v * h_paths_count
    
    # --- Final Calculation ---
    total_paths = total_case1 + total_case2
    
    print("Step-by-step calculation of distinct Prime Paths:")
    print("-" * 50)

    print("Case 1: 1 Horizontal Move & 3 Vertical Moves")
    print(f"Number of 3-move vertical sub-paths (1 -> y1 -> y2 -> 7): {v_paths_count}")
    print(f"Possible placements for the horizontal move: {num_placements_for_h}")
    print(f"Subtotal for Case 1 = {num_placements_for_h} * {v_paths_count} = {total_case1}")
    
    print("\nCase 2: 3 Horizontal Moves & 1 Vertical Move")
    print(f"Number of 3-move horizontal sub-paths (1 -> x1 -> x2 -> 5): {h_paths_count}")
    print(f"Possible placements for the vertical move: {num_placements_for_v}")
    print(f"Subtotal for Case 2 = {num_placements_for_v} * {h_paths_count} = {total_case2}")

    print("-" * 50)
    print("Total Number of Distinct Paths:")
    print(f"{total_case1} (from Case 1) + {total_case2} (from Case 2) = {total_paths}")


solve_prime_path()
<<<104>>>