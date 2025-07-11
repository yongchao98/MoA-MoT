import numpy as np
import math
from collections import Counter

def solve_scanner_placement():
    """
    Solves the museum scanner placement problem using a greedy algorithm.
    """
    # --- Configuration ---
    ROOM_W = 140
    ROOM_H = 110
    COVERAGE_TARGET = 0.88
    
    SCANNERS = {
        'C2': {'shape': 'circle', 'range': 20, 'cost': 20000}, # radius 20m
        'C1': {'shape': 'circle', 'range': 5,  'cost': 1600},  # diameter 10m -> radius 5m
        'R1': {'shape': 'square', 'range': 10, 'cost': 2000},  # side 10m
    }

    # --- Grid Setup ---
    # We use a 1x1m grid. Dimensions are +1 to include the boundaries 0 to 140/110.
    grid_h = ROOM_H + 1
    grid_w = ROOM_W + 1
    total_points = grid_w * grid_h
    target_points_to_cover = math.ceil(total_points * COVERAGE_TARGET)

    # Potential scanner center locations (multiples of 5m)
    locations_x = range(0, ROOM_W + 1, 5)
    locations_y = range(0, ROOM_H + 1, 5)
    potential_locations = [(x, y) for x in locations_x for y in locations_y]

    print("Step 1: Pre-computing coverage masks for all possible scanner placements...")
    # Use coordinate matrices for efficient, vectorized calculations
    yy, xx = np.ogrid[:grid_h, :grid_w]
    
    precomputed_masks = {}
    for name, props in SCANNERS.items():
        shape = props['shape']
        s_range = props['range']
        for loc in potential_locations:
            cx, cy = loc
            if shape == 'circle':
                mask = (xx - cx)**2 + (yy - cy)**2 <= s_range**2
            elif shape == 'square':
                half_side = s_range / 2
                mask = np.logical_and(np.abs(xx - cx) <= half_side, np.abs(yy - cy) <= half_side)
            precomputed_masks[(name, loc)] = mask
    print("Pre-computation complete.")

    # --- Greedy Algorithm ---
    print("\nStep 2: Running greedy algorithm to find the best scanner placement...")
    coverage_grid = np.full((grid_h, grid_w), False)
    placed_scanners = []
    total_cost = 0
    total_covered_points = 0

    while total_covered_points < target_points_to_cover:
        best_choice = None
        best_value = -1  # "value" is new_points / cost
        
        # Find the best scanner to add in this iteration
        for (name, loc), mask in precomputed_masks.items():
            # Skip if this scanner placement doesn't offer any new coverage
            if np.all(coverage_grid[mask]):
                continue

            newly_covered_mask = np.logical_and(mask, ~coverage_grid)
            num_new_points = np.sum(newly_covered_mask)
            
            if num_new_points == 0:
                continue
                
            cost = SCANNERS[name]['cost']
            value = num_new_points / cost
            
            if value > best_value:
                best_value = value
                best_choice = {'name': name, 'location': loc, 'cost': cost}
                
        if best_choice is None:
            print("No further improvement possible. Stopping.")
            break
            
        # Add the best scanner found to our solution
        name = best_choice['name']
        loc = best_choice['location']
        cost = best_choice['cost']
        
        mask_to_add = precomputed_masks[(name, loc)]
        coverage_grid = np.logical_or(coverage_grid, mask_to_add)
        
        placed_scanners.append(best_choice)
        total_cost += cost
        total_covered_points = np.sum(coverage_grid)

    print("Algorithm finished.")
    
    # --- Final Output ---
    print("\n--- Optimization Results ---")
    final_coverage_perc = total_covered_points / total_points
    
    print(f"Target Coverage: {COVERAGE_TARGET:.2%}")
    print(f"Achieved Coverage: {final_coverage_perc:.2%} ({total_covered_points}/{total_points} points)")
    print(f"Total Scanners Used: {len(placed_scanners)}")

    scanner_counts = Counter(s['name'] for s in placed_scanners)
    equation_parts = []
    for name, count in sorted(scanner_counts.items()):
        cost = SCANNERS[name]['cost']
        equation_parts.append(f"({count} * {cost})")
    
    equation_str = f"{' + '.join(equation_parts)} = {total_cost}"
    print(f"\nFinal Cost Equation: {equation_str}")
    
    print(f"\nThe optimal total cost is {total_cost}.")
    print(f"<<<{total_cost}>>>")


if __name__ == '__main__':
    solve_scanner_placement()
