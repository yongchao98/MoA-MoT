import numpy as np
import math

def solve_scanner_placement():
    """
    Solves the scanner placement optimization problem using a greedy algorithm.
    """
    # 1. Define Environment and Scanners
    ROOM_WIDTH = 140
    ROOM_HEIGHT = 110
    GRID_WIDTH = ROOM_WIDTH + 1
    GRID_HEIGHT = ROOM_HEIGHT + 1
    TOTAL_POINTS = GRID_WIDTH * GRID_HEIGHT

    TARGET_COVERAGE_RATIO = 0.88
    TARGET_COVERAGE_POINTS = TOTAL_POINTS * TARGET_COVERAGE_RATIO

    scanners = {
        'C2': {'radius': 20, 'side': None, 'cost': 20000, 'shape': 'circle'},
        'C1': {'radius': 5,  'side': None, 'cost': 1600,  'shape': 'circle'},
        'R1': {'side': 10, 'radius': None, 'cost': 2000,  'shape': 'square'}
    }

    # 2. Setup the grid and possible placements
    center_x_coords = range(0, GRID_WIDTH, 5)
    center_y_coords = range(0, GRID_HEIGHT, 5)
    possible_centers = [(cx, cy) for cx in center_x_coords for cy in center_y_coords]

    # Use a (height, width) or (rows, cols) convention for the numpy array
    coverage_grid = np.zeros((GRID_HEIGHT, GRID_WIDTH), dtype=bool)
    
    # Pre-create coordinate grids for fast calculations
    Y_coords, X_coords = np.mgrid[0:GRID_HEIGHT, 0:GRID_WIDTH]

    # 3. Main Greedy Loop
    total_cost = 0
    placed_scanners = []
    current_coverage_points = 0
    iteration = 0

    print("Starting optimization process...\n")
    
    while current_coverage_points < TARGET_COVERAGE_POINTS:
        iteration += 1
        best_option = {
            'effectiveness': -1,
            'name': None,
            'center': None,
            'cost': -1,
            'new_coverage_mask': None,
            'new_points_count': 0
        }

        # Find the most cost-effective scanner to place next
        for name, props in scanners.items():
            for center in possible_centers:
                cx, cy = center
                
                # Create a boolean mask representing the scanner's coverage area
                if props['shape'] == 'circle':
                    r_sq = props['radius']**2
                    mask = (X_coords - cx)**2 + (Y_coords - cy)**2 <= r_sq
                else: # square
                    half_side = props['side'] / 2
                    mask = (np.abs(X_coords - cx) <= half_side) & (np.abs(Y_coords - cy) <= half_side)

                # Find points that this scanner would cover which are not already covered
                new_coverage_mask = mask & ~coverage_grid
                new_points_count = np.sum(new_coverage_mask)

                if new_points_count == 0:
                    continue

                effectiveness = new_points_count / props['cost']

                if effectiveness > best_option['effectiveness']:
                    best_option.update({
                        'effectiveness': effectiveness,
                        'name': name,
                        'center': center,
                        'cost': props['cost'],
                        'new_coverage_mask': new_coverage_mask,
                        'new_points_count': new_points_count
                    })
        
        if best_option['name'] is None:
            print("Warning: No further placements can improve coverage. Exiting.")
            break
        
        # Place the best scanner from this iteration
        scanner_name = best_option['name']
        scanner_center = best_option['center']
        scanner_cost = best_option['cost']
        
        placed_scanners.append({'name': scanner_name, 'center': scanner_center, 'cost': scanner_cost})
        total_cost += scanner_cost
        
        coverage_grid |= best_option['new_coverage_mask']
        current_coverage_points = np.sum(coverage_grid)
        
        # Uncomment the following line to see step-by-step progress
        # print(f"Iter {iteration}: Placing {scanner_name} at {scanner_center}. New Coverage: {best_option['new_points_count']} pts. Coverage: {current_coverage_points/TOTAL_POINTS:.2%}")

    # 4. Final Output
    print("\n--- Optimization Complete ---\n")

    final_counts = {'C2': 0, 'C1': 0, 'R1': 0}
    for scanner in placed_scanners:
        final_counts[scanner['name']] += 1

    nc2 = final_counts['C2']
    nc1 = final_counts['C1']
    nr1 = final_counts['R1']
    cc2 = scanners['C2']['cost']
    cc1 = scanners['C1']['cost']
    cr1 = scanners['R1']['cost']
    final_coverage_ratio = current_coverage_points / TOTAL_POINTS

    print(f"Final Coverage Achieved: {final_coverage_ratio:.2%} (Target: {TARGET_COVERAGE_RATIO:.2%})")
    print("\nOptimal Scanner Configuration Found:")
    print(f"- C2 Scanners ({cc2} cost): {nc2}")
    print(f"- C1 Scanners ({cc1} cost): {nc1}")
    print(f"- R1 Scanners ({cr1} cost): {nr1}")
    
    print("\nFinal Equation for Total Cost:")
    print(f"{nc2} * {cc2} + {nc1} * {cc1} + {nr1} * {cr1} = {total_cost}")
    
    print(f"\nOptimal Total Cost: {total_cost}")
    
    # Required final answer format
    print(f"\n<<<{total_cost}>>>")

if __name__ == '__main__':
    solve_scanner_placement()