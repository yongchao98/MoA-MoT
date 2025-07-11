import numpy as np
import math

def solve_scanner_placement():
    """
    Solves the museum scanner placement problem using a greedy algorithm.
    """
    # --- Problem Definition ---
    ROOM_WIDTH = 140
    ROOM_HEIGHT = 110
    COVERAGE_TARGET = 0.88

    SCANNERS = {
        'C2': {'shape': 'circle', 'radius': 20, 'cost': 20000},
        'C1': {'shape': 'circle', 'radius': 5,  'cost': 1600}, # 10m diameter
        'R1': {'shape': 'square', 'side': 10,   'cost': 2000}
    }

    # --- Discretization Setup ---
    # Grid for checking coverage (1m resolution)
    GRID_STEP = 1
    grid_x = np.arange(0, ROOM_WIDTH + GRID_STEP, GRID_STEP)
    grid_y = np.arange(0, ROOM_HEIGHT + GRID_STEP, GRID_STEP)
    GRID_COORDS_X, GRID_COORDS_Y = np.meshgrid(grid_x, grid_y, indexing='ij')

    TOTAL_POINTS = len(grid_x) * len(grid_y)
    TARGET_COVERED_POINTS = math.ceil(TOTAL_POINTS * COVERAGE_TARGET)

    # Grid for placing scanners (5m resolution)
    PLACEMENT_STEP = 5
    placement_x = np.arange(0, ROOM_WIDTH + PLACEMENT_STEP, PLACEMENT_STEP)
    placement_y = np.arange(0, ROOM_HEIGHT + PLACEMENT_STEP, PLACEMENT_STEP)

    print("Starting optimization...")
    print(f"Room Area: {ROOM_WIDTH}x{ROOM_HEIGHT} = {ROOM_WIDTH * ROOM_HEIGHT} m^2")
    print(f"Target Coverage: {COVERAGE_TARGET:.2f} ({TARGET_COVERED_POINTS} of {TOTAL_POINTS} grid points)")
    print("-" * 40)
    
    # --- Memoization cache for scanner masks ---
    mask_cache = {}
    def get_coverage_mask(scanner_name, props, center_x, center_y):
        cache_key = (scanner_name, center_x, center_y)
        if cache_key in mask_cache:
            return mask_cache[cache_key]

        if props['shape'] == 'circle':
            mask = (GRID_COORDS_X - center_x)**2 + (GRID_COORDS_Y - center_y)**2 <= props['radius']**2
        elif props['shape'] == 'square':
            half_side = props['side'] / 2
            mask = (np.abs(GRID_COORDS_X - center_x) <= half_side) & \
                   (np.abs(GRID_COORDS_Y - center_y) <= half_side)
        
        mask_cache[cache_key] = mask
        return mask

    # --- Main Greedy Algorithm ---
    placed_scanners = []
    total_cost = 0
    covered_grid = np.zeros((len(grid_x), len(grid_y)), dtype=bool)
    covered_points_count = 0

    while covered_points_count < TARGET_COVERED_POINTS:
        best_candidate = None
        max_effectiveness = -1

        # Find the most cost-effective scanner to add
        for p_x in placement_x:
            for p_y in placement_y:
                for name, props in SCANNERS.items():
                    scanner_mask = get_coverage_mask(name, props, p_x, p_y)
                    newly_covered_mask = np.logical_and(scanner_mask, np.logical_not(covered_grid))
                    newly_covered_count = np.sum(newly_covered_mask)

                    if newly_covered_count > 0:
                        effectiveness = newly_covered_count / props['cost']
                        if effectiveness > max_effectiveness:
                            max_effectiveness = effectiveness
                            best_candidate = {
                                'name': name,
                                'cost': props['cost'],
                                'new_mask': newly_covered_mask
                            }
        
        if best_candidate is None:
            print("Warning: No more coverage can be added. Halting.")
            break

        # Add the best scanner found
        placed_scanners.append(best_candidate)
        total_cost += best_candidate['cost']
        covered_grid = np.logical_or(covered_grid, best_candidate['new_mask'])
        covered_points_count = np.sum(covered_grid)

    # --- Final Output ---
    print("--- Optimization Complete ---")
    
    scanner_counts = {name: 0 for name in SCANNERS}
    for scanner in placed_scanners:
        scanner_counts[scanner['name']] += 1
        
    final_coverage = covered_points_count / TOTAL_POINTS
    
    print(f"\nFinal Coverage Achieved: {final_coverage:.4f}")
    print(f"Total Scanners Placed: {len(placed_scanners)}")
    for name, count in scanner_counts.items():
        if count > 0:
            print(f"  - Type {name}: {count}")

    equation_parts = []
    for name, count in sorted(scanner_counts.items(), key=lambda item: SCANNERS[item[0]]['cost'], reverse=True):
        if count > 0:
            cost = SCANNERS[name]['cost']
            equation_parts.append(f"{count} * {cost}")
    
    equation_str = " + ".join(equation_parts)
    
    print("\nFinal Cost Calculation:")
    print(f"{equation_str} = {total_cost}")
    
    return total_cost

if __name__ == '__main__':
    final_cost = solve_scanner_placement()
    print(f"\n<<<The optimal total cost is: {final_cost}>>>")