import numpy as np

def solve_scanner_placement():
    """
    Solves the museum scanner placement problem using a greedy algorithm.
    """
    # --- Problem Definition ---
    ROOM_W = 140  # meters
    ROOM_H = 110  # meters
    GRID_STEP = 5  # meters for scanner placement
    COVERAGE_TARGET_RATIO = 0.88

    SCANNERS = {
        'C2': {'shape': 'circle', 'radius': 20, 'cost': 20000},
        'C1': {'shape': 'circle', 'radius': 5, 'cost': 1600}, # 10m diameter
        'R1': {'shape': 'square', 'side': 10, 'cost': 2000}
    }

    # --- Grid and Area Calculation ---
    # Use a 1x1 meter resolution for the coverage map
    map_h = ROOM_H + 1
    map_w = ROOM_W + 1
    total_area_points = map_w * map_h
    target_coverage_points = int(np.ceil(total_area_points * COVERAGE_TARGET_RATIO))

    # --- Algorithm Implementation ---

    # 1. Prepare for calculation
    # Create the coordinate grid for the entire room
    gy, gx = np.ogrid[0:map_h, 0:map_w]

    # Get all possible scanner placement center coordinates
    placement_coords = []
    for x in range(0, ROOM_W + 1, GRID_STEP):
        for y in range(0, ROOM_H + 1, GRID_STEP):
            placement_coords.append((x, y))

    # 2. Pre-calculate all possible scanner coverage masks to speed up the main loop
    all_possible_scanners = []
    for name, properties in SCANNERS.items():
        cost = properties['cost']
        for cx, cy in placement_coords:
            if properties['shape'] == 'circle':
                radius = properties['radius']
                mask = (gx - cx)**2 + (gy - cy)**2 <= radius**2
            elif properties['shape'] == 'square':
                half_side = properties['side'] / 2
                mask = (np.abs(gx - cx) <= half_side) & (np.abs(gy - cy) <= half_side)

            all_possible_scanners.append({
                'name': name,
                'center': (cx, cy),
                'cost': cost,
                'mask': mask,
                'area': np.sum(mask) # Store the total area for reference
            })

    # 3. Run the Greedy Algorithm
    coverage_map = np.zeros((map_h, map_w), dtype=bool)
    total_cost = 0
    chosen_scanners = []

    print(f"Goal: Cover at least {target_coverage_points} points ({COVERAGE_TARGET_RATIO:%}) of the {total_area_points} total points.")
    print("Starting greedy optimization...\n")

    while np.sum(coverage_map) < target_coverage_points:
        best_scanner_info = None
        max_efficiency = -1
        
        # This map shows which points are not yet covered
        uncovered_map = ~coverage_map

        # Find the most cost-effective scanner in this iteration
        for scanner in all_possible_scanners:
            # Skip if scanner provides 0 new coverage
            if not np.any(scanner['mask'] & uncovered_map):
                continue

            # Calculate new coverage provided by this scanner
            new_points = np.sum(scanner['mask'] & uncovered_map)
            
            efficiency = new_points / scanner['cost']
            
            if efficiency > max_efficiency:
                max_efficiency = efficiency
                best_scanner_info = scanner

        if best_scanner_info is None:
            print("No more coverage can be added, but target not reached.")
            break

        # Add the best scanner to our solution
        chosen_scanners.append(best_scanner_info)
        coverage_map |= best_scanner_info['mask'] # Update the coverage map
        total_cost += best_scanner_info['cost']

    # 4. Final Output
    print("--- Optimization Complete ---\n")
    
    # Tally the results
    scanner_counts = {'C2': 0, 'C1': 0, 'R1': 0}
    for scanner in chosen_scanners:
        scanner_counts[scanner['name']] += 1

    final_coverage_points = np.sum(coverage_map)
    final_coverage_ratio = final_coverage_points / total_area_points

    print(f"Final Coverage: {final_coverage_points} points ({final_coverage_ratio:.2%})")
    print(f"Final Cost: {total_cost}\n")

    # Print the final equation
    print("Final Cost Breakdown:")
    parts = []
    if scanner_counts['C2'] > 0:
        parts.append(f"{scanner_counts['C2']} * {SCANNERS['C2']['cost']}")
    if scanner_counts['C1'] > 0:
        parts.append(f"{scanner_counts['C1']} * {SCANNERS['C1']['cost']}")
    if scanner_counts['R1'] > 0:
        parts.append(f"{scanner_counts['R1']} * {SCANNERS['R1']['cost']}")
    
    final_equation = " + ".join(parts) + f" = {total_cost}"
    print(final_equation)
    
    return total_cost

if __name__ == '__main__':
    optimal_cost = solve_scanner_placement()
    # The final answer format is handled outside the code block per instructions
    # print(f"\n<<<${optimal_cost}>>>")

# To run, execute the function
solve_scanner_placement()