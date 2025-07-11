import numpy as np

def solve_scanner_placement():
    """
    Uses a greedy algorithm to find a cost-effective scanner arrangement
    to cover a target percentage of a room.
    """
    # 1. Define constants and setup
    ROOM_W, ROOM_H = 140, 110
    TARGET_COVERAGE_RATIO = 0.88
    TOTAL_GRID_POINTS = ROOM_W * ROOM_H
    TARGET_POINTS = int(TOTAL_GRID_POINTS * TARGET_COVERAGE_RATIO)

    SCANNERS_DEF = {
        'C2': {'shape': 'circle', 'param': 20, 'cost': 20000}, # param is radius
        'C1': {'shape': 'circle', 'param': 5,  'cost': 1600},  # param is radius (from 10m diameter)
        'R1': {'shape': 'square', 'param': 10, 'cost': 2000}   # param is side
    }

    print("Problem Setup:")
    print(f"Room Area: {TOTAL_GRID_POINTS} m^2")
    print(f"Target Coverage: {TARGET_COVERAGE_RATIO:.0%} ({TARGET_POINTS} m^2)")
    print("-" * 30)

    # 2. Generate all possible placement coordinates and scanner options
    placement_coords = []
    for x in range(0, ROOM_W + 1, 5):
        for y in range(0, ROOM_H + 1, 5):
            placement_coords.append((x, y))

    # Grid for coverage calculation (1x1m resolution)
    yy, xx = np.mgrid[0:ROOM_H, 0:ROOM_W]

    # 3. Pre-compute coverage masks for every possible scanner placement
    potential_scanners = []
    for name, spec in SCANNERS_DEF.items():
        for center_x, center_y in placement_coords:
            if spec['shape'] == 'circle':
                radius = spec['param']
                # Check if grid point's integer coordinates are within the circle
                mask = (xx - center_x)**2 + (yy - center_y)**2 <= radius**2
            elif spec['shape'] == 'square':
                side = spec['param']
                half_side = side / 2
                x_min, x_max = center_x - half_side, center_x + half_side
                y_min, y_max = center_y - half_side, center_y + half_side
                # Check if grid point's integer coordinates are within the square
                mask = (xx >= x_min) & (xx < x_max) & (yy >= y_min) & (yy < y_max)
            
            potential_scanners.append({
                'name': name,
                'center': (center_x, center_y),
                'cost': spec['cost'],
                'coverage_mask': mask
            })

    # 4. Run the greedy algorithm
    main_coverage_grid = np.zeros((ROOM_H, ROOM_W), dtype=bool)
    total_cost = 0
    chosen_scanners = []
    
    # Use a set of available indices for efficient removal
    available_indices = set(range(len(potential_scanners)))
    used_locations = set()

    while np.sum(main_coverage_grid) < TARGET_POINTS:
        best_effectiveness = -1
        best_scanner_info = None
        best_scanner_idx = -1
        
        # In each iteration, find the most cost-effective scanner to add
        for i in available_indices:
            scanner = potential_scanners[i]
            
            # Skip if a scanner has already been placed at this location
            if scanner['center'] in used_locations:
                continue

            # Calculate how many *new* points this scanner covers
            # Use bitwise operators for speed
            newly_covered_mask = scanner['coverage_mask'] & ~main_coverage_grid
            new_coverage_count = np.sum(newly_covered_mask)

            if new_coverage_count > 0:
                effectiveness = new_coverage_count / scanner['cost']
                if effectiveness > best_effectiveness:
                    best_effectiveness = effectiveness
                    best_scanner_info = scanner
                    best_scanner_idx = i

        if best_scanner_info is None:
            print("Warning: No more scanners can provide additional coverage. Stopping.")
            break
        
        # Add the best scanner found to the solution
        chosen_scanners.append(best_scanner_info)
        total_cost += best_scanner_info['cost']
        main_coverage_grid |= best_scanner_info['coverage_mask']
        
        # Mark this location as used so no other scanner can be placed there
        chosen_location = best_scanner_info['center']
        used_locations.add(chosen_location)

    # 5. Print the final results
    print("Solution Found:")
    final_coverage_pts = np.sum(main_coverage_grid)
    final_coverage_ratio = final_coverage_pts / TOTAL_GRID_POINTS
    print(f"Achieved Coverage: {final_coverage_pts} m^2 ({final_coverage_ratio:.2%})")
    print(f"Total Scanners Used: {len(chosen_scanners)}")
    
    # Print the equation as requested
    costs = [s['cost'] for s in chosen_scanners]
    cost_equation = " + ".join(map(str, costs))
    print("\nCost Equation:")
    print(f"{cost_equation} = {total_cost}")

    print("\nOptimal Total Cost:")
    print(total_cost)
    return total_cost

# Execute the function and capture the final cost for the <<<answer>>> format
final_cost = solve_scanner_placement()
print(f"\n<<<${final_cost}>>>")