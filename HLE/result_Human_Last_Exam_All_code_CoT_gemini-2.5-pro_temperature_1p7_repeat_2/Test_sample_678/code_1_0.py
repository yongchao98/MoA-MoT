import numpy as np
import sys
from collections import Counter

def solve_scanner_placement():
    """
    Uses a greedy algorithm to find a cost-effective scanner placement to achieve target coverage.
    """
    # --- Problem Setup ---
    ROOM_W, ROOM_H = 140, 110
    GRID_STEP = 5
    TARGET_COVERAGE_RATIO = 0.88

    TOTAL_AREA = ROOM_W * ROOM_H
    TARGET_AREA = TOTAL_AREA * TARGET_COVERAGE_RATIO

    SCANNERS = {
        'C2': {'shape': 'circle', 'param': 20, 'cost': 20000},  # param is radius
        'C1': {'shape': 'circle', 'param': 5,  'cost': 1600},   # param is radius (10m diameter)
        'R1': {'shape': 'square', 'param': 10, 'cost': 2000}   # param is side length
    }

    # --- Pre-computation ---
    # Create a grid of coordinates to avoid repeated calculations
    y_coords, x_coords = np.ogrid[0:ROOM_H, 0:ROOM_W]

    # Generate all possible scanner placement locations
    possible_x_centers = range(0, ROOM_W + 1, GRID_STEP)
    possible_y_centers = range(0, ROOM_H + 1, GRID_STEP)
    
    # Pre-calculate the coverage mask for every possible scanner and location
    print("Pre-calculating all possible scanner coverage masks...")
    scanner_masks = {}
    for s_type, s_info in SCANNERS.items():
        cost = s_info['cost']
        shape = s_info['shape']
        param = s_info['param']
        
        for cx in possible_x_centers:
            for cy in possible_y_centers:
                if shape == 'circle':
                    mask = (x_coords - cx)**2 + (y_coords - cy)**2 <= param**2
                elif shape == 'square':
                    half_side = param / 2.0
                    mask = (x_coords >= cx - half_side) & (x_coords < cx + half_side) & \
                           (y_coords >= cy - half_side) & (y_coords < cy + half_side)
                
                # We only need to store masks that cover at least one pixel
                if np.any(mask):
                    scanner_masks[(s_type, cx, cy)] = (mask, cost)
    print("Pre-calculation finished.\n")

    # --- Greedy Algorithm Execution ---
    coverage_grid = np.zeros((ROOM_H, ROOM_W), dtype=bool)
    total_cost = 0
    placements = []
    
    iteration = 1
    while np.sum(coverage_grid) < TARGET_AREA:
        uncovered_grid = ~coverage_grid
        
        best_option = None
        max_effectiveness = -1

        for (s_type, cx, cy), (mask, cost) in scanner_masks.items():
            # Calculate the new area this scanner would cover
            newly_covered_area = np.sum(mask & uncovered_grid)

            if newly_covered_area == 0:
                continue

            effectiveness = newly_covered_area / cost
            
            if effectiveness > max_effectiveness:
                max_effectiveness = effectiveness
                best_option = {
                    'type': s_type,
                    'center': (cx, cy),
                    'cost': cost,
                    'new_area': newly_covered_area,
                    'mask': mask
                }

        if best_option is None:
            print("No more effective placements can be found. Stopping.")
            break
        
        # Add the best scanner from this iteration
        placements.append(best_option)
        total_cost += best_option['cost']
        coverage_grid |= best_option['mask']
        current_covered_area = np.sum(coverage_grid)
        
        print(f"Step {iteration}: Placed {best_option['type']} at {best_option['center']}. "
              f"Coverage: {current_covered_area / TOTAL_AREA:.2%}. Cost: {total_cost}")
        iteration += 1
    
    # --- Final Results ---
    print("\n--- Optimization Complete ---")
    final_coverage = np.sum(coverage_grid)
    final_coverage_ratio = final_coverage / TOTAL_AREA
    
    print(f"Final Coverage: {final_coverage} / {TOTAL_AREA} m^2 ({final_coverage_ratio:.4f})")
    print(f"Target Coverage: {TARGET_AREA} / {TOTAL_AREA} m^2 ({TARGET_COVERAGE_RATIO:.4f})")
    
    scanner_counts = Counter(p['type'] for p in placements)
    n_c2 = scanner_counts['C2']
    n_c1 = scanner_counts['C1']
    n_r1 = scanner_counts['R1']
    
    print("\nScanner Count:")
    print(f" - C2 (20m radius): {n_c2}")
    print(f" - C1 (10m diameter): {n_c1}")
    print(f" - R1 (10m square): {n_r1}")
    
    cost_c2 = SCANNERS['C2']['cost']
    cost_c1 = SCANNERS['C1']['cost']
    cost_r1 = SCANNERS['R1']['cost']

    equation = (
        f"Total Cost = ({n_c2} * {cost_c2}) + ({n_c1} * {cost_c1}) + ({n_r1} * {cost_r1}) = {total_cost}"
    )
    print("\nFinal Cost Calculation:")
    print(equation)

    print("\n<<<" + str(total_cost) + ">>>")


if __name__ == '__main__':
    solve_scanner_placement()