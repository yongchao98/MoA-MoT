import numpy as np
import math

def solve_scanner_placement():
    """
    Solves the museum scanner placement problem using a greedy algorithm.
    """
    # 1. Define constants and setup
    ROOM_W = 140  # Room width in meters
    ROOM_H = 110  # Room height in meters
    TARGET_COVERAGE_RATIO = 0.88

    SCANNERS = {
        'C2': {'type': 'circle', 'radius': 20, 'cost': 20000},
        'C1': {'type': 'circle', 'radius': 5,  'cost': 1600}, # 10m diameter
        'R1': {'type': 'square', 'side': 10, 'cost': 2000}
    }

    # 2. Setup the grid and coverage targets
    total_grid_points = ROOM_W * ROOM_H
    target_covered_points = math.ceil(total_grid_points * TARGET_COVERAGE_RATIO)

    # A boolean grid representing the room. True means covered.
    # We use (y, x) indexing to align with numpy's (row, col) convention.
    covered_grid = np.zeros((ROOM_H, ROOM_W), dtype=bool)

    # Define the discrete set of possible center locations for scanners
    possible_cx = range(0, ROOM_W + 1, 5)
    possible_cy = range(0, ROOM_H + 1, 5)
    
    # 3. Main greedy loop
    total_cost = 0
    placed_scanners = []
    
    print("Starting optimization...")
    print(f"Room Area: {total_grid_points} sq.m.")
    print(f"Target Coverage: {TARGET_COVERAGE_RATIO * 100}% ({target_covered_points} sq.m.)\n")

    iteration = 0
    while np.sum(covered_grid) < target_covered_points:
        iteration += 1
        best_placement = None
        best_marginal_gain = -1.0

        # Iterate through every possible scanner at every possible location
        for name, info in SCANNERS.items():
            cost = info['cost']
            for cx in possible_cx:
                for cy in possible_cy:
                    
                    newly_covered_coords = []
                    # Determine the coordinates this scanner would cover
                    if info['type'] == 'circle':
                        r = info['radius']
                        r_squared = r**2
                        # Check points within the scanner's bounding box for efficiency
                        y_min, y_max = max(0, cy - r), min(ROOM_H, cy + r + 1)
                        x_min, x_max = max(0, cx - r), min(ROOM_W, cx + r + 1)
                        
                        for y in range(y_min, y_max):
                            for x in range(x_min, x_max):
                                # If point is not already covered and is inside the circle
                                if not covered_grid[y, x] and (x - cx)**2 + (y - cy)**2 <= r_squared:
                                    newly_covered_coords.append((y, x))
                    
                    elif info['type'] == 'square':
                        s_half = info['side'] / 2.0
                        # Get integer range of covered coordinates
                        y_min, y_max = math.ceil(cy - s_half), math.floor(cy + s_half)
                        x_min, x_max = math.ceil(cx - s_half), math.floor(cx + s_half)
                        
                        # Clamp to room boundaries
                        y_min, y_max = max(0, y_min), min(ROOM_H, y_max)
                        x_min, x_max = max(0, x_min), min(ROOM_W, x_max)
                        
                        for y in range(y_min, y_max):
                            for x in range(x_min, x_max):
                                # If point is not already covered
                                if not covered_grid[y, x]:
                                    newly_covered_coords.append((y, x))
                    
                    # Calculate marginal gain for this potential placement
                    new_points_count = len(newly_covered_coords)
                    if new_points_count > 0:
                        marginal_gain = new_points_count / cost
                        if marginal_gain > best_marginal_gain:
                            best_marginal_gain = marginal_gain
                            best_placement = {
                                'name': name,
                                'center': (cx, cy),
                                'cost': cost,
                                'coords': newly_covered_coords,
                                'new_points': new_points_count
                            }

        # If no placement can add new coverage, break the loop
        if best_placement is None:
            print("Warning: Could not find any more placements to improve coverage.")
            break

        # 4. "Place" the best scanner found in this iteration
        name = best_placement['name']
        center = best_placement['center']
        cost = best_placement['cost']
        
        placed_scanners.append({'name': name, 'center': center, 'cost': cost})
        total_cost += cost
        
        # Update the grid to mark the newly covered points
        coords_to_update = best_placement['coords']
        rows, cols = zip(*coords_to_update)
        covered_grid[rows, cols] = True
        
        current_coverage_points = np.sum(covered_grid)
        current_coverage_ratio = current_coverage_points / total_grid_points
        
        print(f"Iteration {iteration}: Placed {name} at {center}. "
              f"Coverage: {current_coverage_ratio:.2%}. "
              f"Cost: {total_cost}")

    # 5. Print final results
    print("\n--- Optimization Complete ---")
    
    n_c2 = sum(1 for s in placed_scanners if s['name'] == 'C2')
    n_c1 = sum(1 for s in placed_scanners if s['name'] == 'C1')
    n_r1 = sum(1 for s in placed_scanners if s['name'] == 'R1')
    
    cost_c2 = SCANNERS['C2']['cost']
    cost_c1 = SCANNERS['C1']['cost']
    cost_r1 = SCANNERS['R1']['cost']

    print(f"Total Scanners Placed: {len(placed_scanners)}")
    print(f" - C2 Scanners (r=20m): {n_c2}")
    print(f" - C1 Scanners (r=5m): {n_c1}")
    print(f" - R1 Scanners (s=10m): {n_r1}")
    
    final_coverage_points = np.sum(covered_grid)
    final_coverage_ratio = final_coverage_points / total_grid_points
    print(f"\nFinal Coverage: {final_coverage_points} / {total_grid_points} sq.m. = {final_coverage_ratio:.2%}")
    
    print("\nFinal Cost Calculation:")
    print(f"({n_c2} * {cost_c2}) + ({n_c1} * {cost_c1}) + ({n_r1} * {cost_r1}) = {total_cost}")

    print("\nOptimal Total Cost:")
    print(total_cost)
    
    return total_cost

if __name__ == '__main__':
    final_cost = solve_scanner_placement()
    print(f"\n<<<{final_cost}>>>")