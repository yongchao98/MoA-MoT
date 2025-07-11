import numpy as np
import math

def solve_scanner_placement():
    """
    Formulates and solves the museum scanner placement problem using a greedy algorithm.
    """
    # 1. SETUP
    
    # Room and Target Coverage
    ROOM_WIDTH = 140
    ROOM_HEIGHT = 110
    TOTAL_AREA = ROOM_WIDTH * ROOM_HEIGHT
    TARGET_COVERAGE_RATIO = 0.88
    TARGET_COVERAGE_AREA = TARGET_COVERAGE_RATIO * TOTAL_AREA

    # Scanner definitions
    SCANNERS = {
        'C2': {'shape': 'circle', 'radius': 20, 'cost': 20000},
        'C1': {'shape': 'circle', 'radius': 5, 'cost': 1600}, # from 10m diameter
        'R1': {'shape': 'square', 'side': 10, 'cost': 2000}
    }

    # Valid placement points (centers are multiples of 5m)
    placement_xs = range(0, ROOM_WIDTH + 1, 5)
    placement_ys = range(0, ROOM_HEIGHT + 1, 5)
    placement_points = [(x, y) for x in placement_xs for y in placement_ys]

    # Coverage grid (1m x 1m resolution)
    # We model coverage by checking if the center of each 1m^2 cell is covered.
    coverage_grid = np.zeros((ROOM_HEIGHT, ROOM_WIDTH), dtype=bool)
    grid_y_centers, grid_x_centers = np.mgrid[0:ROOM_HEIGHT, 0:ROOM_WIDTH]
    grid_x_centers = grid_x_centers + 0.5
    grid_y_centers = grid_y_centers + 0.5

    # 2. PRE-COMPUTATION of scanner masks

    def create_mask(scanner_type, center_x, center_y):
        """Creates a boolean mask for a given scanner and location."""
        s_info = SCANNERS[scanner_type]
        if s_info['shape'] == 'circle':
            radius = s_info['radius']
            # Check if cell center is within the circle
            dist_sq = (grid_x_centers - center_x)**2 + (grid_y_centers - center_y)**2
            return dist_sq <= radius**2
        elif s_info['shape'] == 'square':
            half_side = s_info['side'] / 2.0
            # Check if cell center is within the square
            return (np.abs(grid_x_centers - center_x) <= half_side) & \
                   (np.abs(grid_y_centers - center_y) <= half_side)

    all_masks = {}
    for name in SCANNERS:
        all_masks[name] = {
            (p_x, p_y): create_mask(name, p_x, p_y)
            for p_x, p_y in placement_points
        }

    # 3. GREEDY ALGORITHM to place scanners

    total_cost = 0
    placed_scanners_count = {'C2': 0, 'C1': 0, 'R1': 0}

    while np.sum(coverage_grid) < TARGET_COVERAGE_AREA:
        best_efficiency = -1
        best_scanner_info = None
        best_mask = None

        # Find the most cost-effective scanner to add in this iteration
        for s_name, s_info in SCANNERS.items():
            cost = s_info['cost']
            for p_coords, mask in all_masks[s_name].items():
                
                # Calculate new coverage (pixels in mask but not yet covered)
                # Using bitwise operations for efficiency
                newly_covered_mask = mask & ~coverage_grid
                newly_covered_area = np.sum(newly_covered_mask)

                if newly_covered_area == 0:
                    continue

                efficiency = newly_covered_area / cost
                if efficiency > best_efficiency:
                    best_efficiency = efficiency
                    best_scanner_info = (s_name, p_coords)
                    best_mask = mask

        if best_scanner_info is None:
            # This case happens if no scanner can add any new coverage
            break

        # Place the best scanner found
        s_name_best, _ = best_scanner_info
        total_cost += SCANNERS[s_name_best]['cost']
        placed_scanners_count[s_name_best] += 1
        
        # Update the master coverage grid
        coverage_grid |= best_mask

    # 4. OUTPUT RESULTS
    
    final_coverage_area = np.sum(coverage_grid)
    final_coverage_ratio = final_coverage_area / TOTAL_AREA

    print("--- Optimization Result ---")
    print(f"Final Coverage: {final_coverage_area} m^2 / {TOTAL_AREA} m^2 ({final_coverage_ratio:.4f})")
    print(f"Target Coverage: >={TARGET_COVERAGE_AREA} m^2 ({TARGET_COVERAGE_RATIO:.2f})")
    print("\nOptimal Scanner Configuration:")
    for name, count in placed_scanners_count.items():
        if count > 0:
            print(f"- Type {name}: {count} units")
    
    print("\nFinal Cost Calculation:")
    cost_parts = []
    for name, count in placed_scanners_count.items():
        if count > 0:
            cost_parts.append(f"{count} * {SCANNERS[name]['cost']}")
    
    equation = f"{' + '.join(cost_parts)} = {total_cost}"
    print(equation)
    
    print(f"\nOptimal Total Cost: {total_cost}")
    print(f"<<<{total_cost}>>>")

# Execute the function
solve_scanner_placement()