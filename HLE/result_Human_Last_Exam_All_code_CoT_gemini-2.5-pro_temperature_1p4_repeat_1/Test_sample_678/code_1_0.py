import numpy as np
import math
from collections import Counter

def solve_scanner_placement():
    """
    Solves the museum scanner placement problem using a greedy algorithm.
    """
    # 1. Define constants and room setup
    ROOM_W, ROOM_H = 140, 110
    TARGET_COVERAGE_RATIO = 0.88
    PLACEMENT_GRID_STEP = 5

    scanners = {
        'C2': {'type': 'circle', 'radius': 20, 'cost': 20000},
        'C1': {'type': 'circle', 'radius': 5,  'cost': 1600}, # 10m diameter
        'R1': {'type': 'square', 'side': 10,   'cost': 2000},
    }

    # Discretize the room into a 1x1m grid
    # Grid dimensions are +1 because we include the boundaries (0 to 140, 0 to 110)
    grid_w, grid_h = ROOM_W + 1, ROOM_H + 1
    total_points = grid_w * grid_h
    target_covered_points = math.floor(total_points * TARGET_COVERAGE_RATIO)

    print(f"Room Area: {ROOM_W}x{ROOM_H}m = {ROOM_W * ROOM_H} m^2")
    print(f"Target Coverage: {TARGET_COVERAGE_RATIO:.0%}")
    print(f"Total Grid Points: {total_points}")
    print(f"Target Covered Points: {target_covered_points}\n")
    
    # 2. Pre-compute potential scanner coverage masks
    
    # Generate possible center locations
    possible_centers = []
    for x in range(0, ROOM_W + 1, PLACEMENT_GRID_STEP):
        for y in range(0, ROOM_H + 1, PLACEMENT_GRID_STEP):
            possible_centers.append((x, y))

    # Create grid coordinates for vector-based calculations
    Y, X = np.ogrid[:grid_h, :grid_w]

    print("Pre-computing scanner coverage masks...")
    precomputed_masks = {}
    for name, spec in scanners.items():
        for cx, cy in possible_centers:
            if spec['type'] == 'circle':
                radius = spec['radius']
                mask = (X - cx)**2 + (Y - cy)**2 <= radius**2
            elif spec['type'] == 'square':
                half_side = spec['side'] / 2
                mask = (np.abs(X - cx) <= half_side) & (np.abs(Y - cy) <= half_side)
            
            # Only store if the mask covers at least one point
            if np.any(mask):
                 precomputed_masks[(name, (cx, cy))] = mask
    print(f"Done. Found {len(precomputed_masks)} potential scanner placements.\n")

    # 3. Run the greedy algorithm
    coverage_map = np.zeros((grid_h, grid_w), dtype=bool)
    total_cost = 0
    placements = []
    
    covered_points_count = 0
    
    while covered_points_count < target_covered_points:
        best_marginal_gain = -1
        best_choice = None
        best_mask = None

        for (name, center), mask in precomputed_masks.items():
            # Calculate newly covered points
            # The '& ~coverage_map' part finds points in the new mask that are not yet covered
            newly_covered_mask = mask & ~coverage_map
            num_new_points = np.sum(newly_covered_mask)

            if num_new_points == 0:
                continue
            
            cost = scanners[name]['cost']
            marginal_gain = num_new_points / cost

            if marginal_gain > best_marginal_gain:
                best_marginal_gain = marginal_gain
                best_choice = (name, center, cost)
                best_mask = newly_covered_mask
        
        if best_choice is None:
            print("No more coverage can be added. Stopping.")
            break
            
        # Add the best scanner from this iteration to the solution
        name, center, cost = best_choice
        placements.append(best_choice)
        total_cost += cost
        coverage_map |= best_mask # Update the master coverage map
        
        new_total_covered = np.sum(coverage_map)
        
        print(f"Added: {name} at {center}. "
              f"Cost: ${cost}. "
              f"New Coverage: {new_total_covered / total_points:.2%}. "
              f"Total Cost: ${total_cost}")
        
        covered_points_count = new_total_covered

    # 4. Final Output
    print("\n--- Optimization Complete ---")

    if not placements:
        print("No scanners were placed.")
        return

    # Summarize the results
    scanner_counts = Counter(p[0] for p in placements)
    
    cost_components = []
    for name, spec in sorted(scanners.items(), key=lambda item: item[1]['cost'], reverse=True):
        count = scanner_counts[name]
        cost = spec['cost']
        print(f"Number of {name} scanners: {count}")
        cost_components.append(f"{count} * {cost}")
        
    final_equation = " + ".join(cost_components)
    print(f"\nFinal Cost Equation: {final_equation} = {total_cost}")

    final_coverage_percent = (covered_points_count / total_points) * 100
    print(f"Final Coverage: {final_coverage_percent:.2f}% ({covered_points_count}/{total_points} points)")
    print(f"Optimal Total Cost: {total_cost}")
    print(f"<<<{total_cost}>>>")

if __name__ == '__main__':
    solve_scanner_placement()