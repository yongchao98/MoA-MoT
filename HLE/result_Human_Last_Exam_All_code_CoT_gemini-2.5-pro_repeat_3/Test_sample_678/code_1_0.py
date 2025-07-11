import numpy as np
import math

def solve_scanner_placement():
    """
    Solves the museum scanner placement problem using a greedy algorithm
    to find a cost-effective solution for covering a specified area.
    """
    # --- Problem Definition ---
    ROOM_W = 140  # meters
    ROOM_H = 110  # meters
    TARGET_COVERAGE_RATIO = 0.88

    SCANNERS = [
        {'name': 'C2', 'cost': 20000, 'type': 'circle', 'param': 20},  # param is radius
        {'name': 'C1', 'cost': 1600,  'type': 'circle', 'param': 5},   # param is radius (10m diameter)
        {'name': 'R1', 'cost': 2000,  'type': 'square', 'param': 10}   # param is side
    ]

    # --- Step 1: Initialize Grid and Targets ---
    print("Step 1: Initializing the room grid and coverage targets...")
    total_area = ROOM_W * ROOM_H
    target_area = total_area * TARGET_COVERAGE_RATIO
    
    # Create a grid representing the room floor. Each cell is 1x1m.
    # coverage_grid[y, x] is True if the cell at (x, y) is covered.
    coverage_grid = np.zeros((ROOM_H, ROOM_W), dtype=bool)
    
    # Create coordinate arrays for vectorized calculations. We use the center of each 1x1m cell.
    grid_y, grid_x = np.mgrid[0.5:ROOM_H, 0.5:ROOM_W]

    # --- Step 2: Define Potential Scanner Locations ---
    print("Step 2: Defining all possible scanner locations (multiples of 5m).")
    potential_x_coords = range(0, ROOM_W + 1, 5)
    potential_y_coords = range(0, ROOM_H + 1, 5)
    potential_locations = [(x, y) for x in potential_x_coords for y in potential_y_coords]

    # --- Step 3: Run the Greedy Algorithm ---
    print("\nStep 3: Running greedy optimization...")
    print("Continuously adding the most cost-effective scanner until coverage target is met.")
    
    total_cost = 0
    current_covered_area = 0
    scanner_counts = {s['name']: 0 for s in SCANNERS}
    
    iteration = 1
    while current_covered_area < target_area:
        best_candidate = None
        min_cost_per_new_area = float('inf')

        # Evaluate every scanner type at every possible location
        for scanner_info in SCANNERS:
            for cx, cy in potential_locations:
                # Calculate the coverage mask for the current candidate scanner
                if scanner_info['type'] == 'circle':
                    mask = (grid_x - cx)**2 + (grid_y - cy)**2 <= scanner_info['param']**2
                else:  # square
                    half_side = scanner_info['param'] / 2.0
                    mask = (np.abs(grid_x - cx) <= half_side) & \
                           (np.abs(grid_y - cy) <= half_side)

                # Identify cells that are covered by this scanner but not yet by others
                newly_covered_mask = mask & ~coverage_grid
                new_area_gain = np.sum(newly_covered_mask)

                if new_area_gain > 0:
                    # Calculate cost-effectiveness
                    cost_effectiveness = scanner_info['cost'] / new_area_gain
                    
                    if cost_effectiveness < min_cost_per_new_area:
                        min_cost_per_new_area = cost_effectiveness
                        best_candidate = {
                            'name': scanner_info['name'],
                            'cost': scanner_info['cost'],
                            'mask': newly_covered_mask,
                            'gain': new_area_gain
                        }

        # After checking all candidates, add the best one to our solution
        if best_candidate:
            scanner_name = best_candidate['name']
            total_cost += best_candidate['cost']
            current_covered_area += best_candidate['gain']
            scanner_counts[scanner_name] += 1
            coverage_grid |= best_candidate['mask']
            
            print(f" Iteration {iteration}: Added {scanner_name} scanner. "
                  f"Total Cost: ${total_cost}, "
                  f"Coverage: {100 * current_covered_area / total_area:.2f}%")
            iteration += 1
        else:
            print("Warning: No more effective scanners found. Stopping.")
            break

    # --- Step 4: Present the Final Answer ---
    print("\nOptimization Complete.")
    print("The proposed solution to achieve at least 88% coverage is:")
    
    num_c2 = scanner_counts['C2']
    num_c1 = scanner_counts['C1']
    num_r1 = scanner_counts['R1']
    
    cost_c2 = next(s['cost'] for s in SCANNERS if s['name'] == 'C2')
    cost_c1 = next(s['cost'] for s in SCANNERS if s['name'] == 'C1')
    cost_r1 = next(s['cost'] for s in SCANNERS if s['name'] == 'R1')

    final_coverage_percent = 100 * current_covered_area / total_area
    print(f"\nFinal Coverage Achieved: {current_covered_area} / {total_area} sq m = {final_coverage_percent:.2f}%")
    
    print("\nOptimal Cost Calculation:")
    print(f"{num_c2} * ${cost_c2} (C2) + {num_c1} * ${cost_c1} (C1) + {num_r1} * ${cost_r1} (R1) = ${total_cost}")

if __name__ == '__main__':
    solve_scanner_placement()