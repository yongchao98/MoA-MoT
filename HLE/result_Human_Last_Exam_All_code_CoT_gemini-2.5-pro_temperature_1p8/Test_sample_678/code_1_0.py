import numpy as np
import math

def solve_scanner_placement():
    """
    Solves the scanner placement problem using a greedy algorithm to minimize cost
    while achieving a target coverage percentage.
    """

    # --- Problem Configuration ---
    ROOM_W_M = 140
    ROOM_H_M = 110
    PLACEMENT_STEP_M = 5
    TARGET_COVERAGE_RATIO = 0.88

    SCANNERS = {
        'C2': {'type': 'circle', 'radius': 20, 'cost': 20000},
        'C1': {'type': 'circle', 'radius': 5,  'cost': 1600},  # 10m diameter
        'R1': {'type': 'square', 'side': 10, 'cost': 2000}
    }

    # --- Grid Setup for Area Calculation (1 cell = 1 m^2) ---
    TOTAL_AREA_M2 = ROOM_W_M * ROOM_H_M
    TARGET_COVERAGE_M2 = TOTAL_AREA_M2 * TARGET_COVERAGE_RATIO

    # Create a meshgrid representing each 1x1m cell in the room.
    # This allows for fast, vectorized calculations of scanner coverage.
    grid_y, grid_x = np.mgrid[0:ROOM_H_M, 0:ROOM_W_M]

    # A cache to store pre-calculated scanner coverage masks to speed up the process.
    mask_cache = {}

    def get_scanner_mask(scanner_name, center_x, center_y):
        """Calculates and returns a boolean mask for a given scanner's coverage."""
        cache_key = (scanner_name, center_x, center_y)
        if cache_key in mask_cache:
            return mask_cache[cache_key]

        scanner_info = SCANNERS[scanner_name]
        mask = np.zeros_like(grid_x, dtype=bool)

        if scanner_info['type'] == 'circle':
            radius = scanner_info['radius']
            # Calculate squared distance from center to every point on the grid
            distance_sq = (grid_x - center_x)**2 + (grid_y - center_y)**2
            mask = distance_sq <= radius**2
        elif scanner_info['type'] == 'square':
            side = scanner_info['side']
            half_side = side / 2
            # Check if grid points fall within the square's boundaries
            mask = (grid_x >= center_x - half_side) & (grid_x < center_x + half_side) & \
                   (grid_y >= center_y - half_side) & (grid_y < center_y + half_side)

        mask_cache[cache_key] = mask
        return mask

    # --- Greedy Algorithm ---
    possible_placements_x = range(0, ROOM_W_M + 1, PLACEMENT_STEP_M)
    possible_placements_y = range(0, ROOM_H_M + 1, PLACEMENT_STEP_M)

    # Main grid to track total area covered so far.
    total_covered_grid = np.zeros_like(grid_x, dtype=bool)
    
    total_cost = 0
    scanner_counts = {name: 0 for name in SCANNERS}

    iteration = 1
    while np.sum(total_covered_grid) < TARGET_COVERAGE_M2:
        best_candidate = {'effectiveness': -1}

        # Evaluate every scanner type at every possible location
        for s_name, s_info in SCANNERS.items():
            for cx in possible_placements_x:
                for cy in possible_placements_y:
                    potential_mask = get_scanner_mask(s_name, cx, cy)
                    
                    # Calculate new area by finding cells covered by this scanner
                    # but not by the already placed ones.
                    newly_covered_mask = potential_mask & ~total_covered_grid
                    newly_covered_area = np.sum(newly_covered_mask)

                    if newly_covered_area == 0:
                        continue  # This option adds no value

                    effectiveness = newly_covered_area / s_info['cost']

                    if effectiveness > best_candidate['effectiveness']:
                        best_candidate = {
                            'effectiveness': effectiveness,
                            'name': s_name,
                            'cost': s_info['cost'],
                            'new_mask': newly_covered_mask
                        }
        
        if best_candidate['effectiveness'] == -1:
            print("Warning: Could not find any more effective placements. Stopping.")
            break

        # Place the best scanner found in this iteration
        total_cost += best_candidate['cost']
        total_covered_grid |= best_candidate['new_mask']
        scanner_counts[best_candidate['name']] += 1
        
        iteration += 1

    # --- Output Results ---
    print("--- Optimization Complete ---")
    final_coverage_m2 = np.sum(total_covered_grid)
    final_coverage_ratio = final_coverage_m2 / TOTAL_AREA_M2

    print(f"Final Coverage: {final_coverage_m2} / {TOTAL_AREA_M2} m^2 = {final_coverage_ratio:.2%}")
    print(f"Target Coverage: {TARGET_COVERAGE_RATIO:.2%}\n")

    cost_components = []
    for name, count in scanner_counts.items():
        if count > 0:
            cost = SCANNERS[name]['cost']
            cost_components.append(f"{count} * {cost}")
    
    equation_str = " + ".join(cost_components)
    
    print(f"Optimal cost calculated from: {equation_str} = {total_cost}")
    print("\nFinal Answer: The optimal total cost is {total_cost}.")
    print(f'<<<{total_cost}>>>')


# Run the optimization
solve_scanner_placement()
