import numpy as np

def solve_scanner_placement():
    """
    Solves the museum scanner placement problem by finding the lowest cost
    configuration that achieves at least 88% coverage.
    """
    # 1. Problem Definition
    ROOM_L = 140
    ROOM_W = 110
    COVERAGE_TARGET_RATIO = 0.88
    PLACEMENT_STEP = 5

    ROOM_AREA = ROOM_L * ROOM_W
    TARGET_COVERAGE_AREA = ROOM_AREA * COVERAGE_TARGET_RATIO

    SCANNERS = {
        'C2': {'shape': 'circle', 'radius': 20, 'cost': 20000},
        'C1': {'shape': 'circle', 'radius': 5,  'cost': 1600},
        'R1': {'shape': 'square', 'side': 10, 'cost': 2000}
    }

    # 2. Setup Grids
    # Create a 1x1m grid to evaluate coverage
    eval_y, eval_x = np.mgrid[0:ROOM_W, 0:ROOM_L]

    # Define all possible placement locations for scanner centers
    placement_xs = range(0, ROOM_L + 1, PLACEMENT_STEP)
    placement_ys = range(0, ROOM_W + 1, PLACEMENT_STEP)

    # 3. Pre-calculate all possible scanner masks to optimize the main loop
    all_possible_scanners = []
    for name, props in SCANNERS.items():
        for cx in placement_xs:
            for cy in placement_ys:
                if props['shape'] == 'circle':
                    mask = (eval_x - cx)**2 + (eval_y - cy)**2 <= props['radius']**2
                else:  # square
                    half_side = props['side'] / 2
                    mask = (np.abs(eval_x - cx) <= half_side) & \
                           (np.abs(eval_y - cy) <= half_side)
                
                # Only consider placements that cover some area within the room
                if np.any(mask):
                    all_possible_scanners.append({
                        'name': name, 'cx': cx, 'cy': cy,
                        'cost': props['cost'], 'mask': mask
                    })

    # 4. Run the Greedy Algorithm
    coverage_map = np.zeros((ROOM_W, ROOM_L), dtype=bool)
    total_cost = 0
    placed_scanners = []
    current_covered_area = 0

    while current_covered_area < TARGET_COVERAGE_AREA:
        best_option = {
            'cost_per_m2': float('inf'),
            'scanner_info': None
        }

        # Find the most cost-effective scanner to add in this iteration
        for scanner in all_possible_scanners:
            # Calculate the new area this scanner would cover
            newly_covered_area = np.sum(scanner['mask'] & ~coverage_map)

            if newly_covered_area > 0:
                cost_per_m2 = scanner['cost'] / newly_covered_area
                if cost_per_m2 < best_option['cost_per_m2']:
                    best_option['cost_per_m2'] = cost_per_m2
                    best_option['scanner_info'] = scanner

        if best_option['scanner_info'] is None:
            print("Warning: No further coverage possible, but target not met.")
            break
        
        # Add the best scanner from this iteration to the solution
        best_scanner = best_option['scanner_info']
        placed_scanners.append(best_scanner)
        total_cost += best_scanner['cost']
        
        # Update the master coverage map
        coverage_map |= best_scanner['mask']
        current_covered_area = np.sum(coverage_map)

        # To speed up, remove the chosen scanner from future consideration
        all_possible_scanners.remove(best_scanner)

    # 5. Print the final results
    print(f"Optimization Complete.")
    print(f"Target Coverage: >={TARGET_COVERAGE_AREA:,.0f} m^2 ({COVERAGE_TARGET_RATIO:.0%})")
    print(f"Achieved Coverage: {current_covered_area:,.0f} m^2 ({current_covered_area / ROOM_AREA:.2%})")
    print(f"Optimal Total Cost: {total_cost}")

    counts = {'C2': 0, 'C1': 0, 'R1': 0}
    for scanner in placed_scanners:
        counts[scanner['name']] += 1

    print("\nFinal Cost Calculation:")
    cost_str = (
        f"{counts['C2']} * {SCANNERS['C2']['cost']} (C2) + "
        f"{counts['C1']} * {SCANNERS['C1']['cost']} (C1) + "
        f"{counts['R1']} * {SCANNERS['R1']['cost']} (R1) = {total_cost}"
    )
    print(cost_str)
    
    return total_cost

if __name__ == '__main__':
    final_cost = solve_scanner_placement()
    print(f"<<<{final_cost}>>>")