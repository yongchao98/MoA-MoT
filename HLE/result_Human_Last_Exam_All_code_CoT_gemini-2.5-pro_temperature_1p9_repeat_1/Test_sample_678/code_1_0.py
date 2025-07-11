import numpy as np
import math

def solve_scanner_placement():
    """
    Solves the museum scanner placement problem using a greedy algorithm
    to find a low-cost solution for covering at least 88% of the area.
    """
    # 1. Constants and Setup
    ROOM_WIDTH = 140  # meters
    ROOM_HEIGHT = 110 # meters
    ROOM_AREA = ROOM_WIDTH * ROOM_HEIGHT
    COVERAGE_GOAL = 0.88
    TARGET_AREA = ROOM_AREA * COVERAGE_GOAL

    # Scanner definitions: (Type, parameter, Cost)
    # For circles, parameter is radius. For squares, parameter is side length.
    SCANNERS = [
        {'type': 'C2', 'param': 20, 'cost': 20000, 'shape': 'circle'},
        {'type': 'C1', 'param': 5,  'cost': 1600,  'shape': 'circle'}, # Diameter 10m -> Radius 5m
        {'type': 'R1', 'param': 10, 'cost': 2000,  'shape': 'square'}
    ]

    # 2. Grid and Placement Points
    placement_x = range(0, ROOM_WIDTH + 1, 5)
    placement_y = range(0, ROOM_HEIGHT + 1, 5)

    # Coordinate grid for area calculation (centers of 1x1m squares)
    yy, xx = np.mgrid[0.5:ROOM_HEIGHT, 0.5:ROOM_WIDTH]

    # 3. Pre-computation of Masks for efficiency
    print("Pre-computing all possible scanner coverage masks...")
    precomputed_masks = {}
    for p_x in placement_x:
        for p_y in placement_y:
            for scanner in SCANNERS:
                key = (scanner['type'], p_x, p_y)
                if scanner['shape'] == 'circle':
                    radius = scanner['param']
                    mask = (xx - p_x)**2 + (yy - p_y)**2 <= radius**2
                elif scanner['shape'] == 'square':
                    side = scanner['param']
                    half_side = side / 2.0
                    mask = (np.abs(xx - p_x) <= half_side) & (np.abs(yy - p_y) <= half_side)
                precomputed_masks[key] = mask
    print("Pre-computation complete.")

    # 4. Greedy Algorithm Loop
    coverage_grid = np.zeros((ROOM_HEIGHT, ROOM_WIDTH), dtype=bool)
    total_cost = 0
    total_area = 0
    scanner_counts = {s['type']: 0 for s in SCANNERS}
    scanner_costs = {s['type']: s['cost'] for s in SCANNERS}

    print("\nStarting optimization...")
    print(f"Room Area: {ROOM_AREA} m^2, Target Coverage: {TARGET_AREA:.2f} m^2 ({COVERAGE_GOAL*100}%)")
    
    step = 1
    while total_area < TARGET_AREA:
        best_effectiveness = -1
        best_choice_key = None
        best_choice_new_area = 0

        # Iterate through all possibilities to find the best next move
        for key, mask in precomputed_masks.items():
            scanner_type = key[0]
            cost = scanner_costs[scanner_type]
            
            # This is faster than creating a new mask for the difference
            new_area = np.sum(mask & ~coverage_grid)

            if new_area > 0:
                effectiveness = new_area / cost
                if effectiveness > best_effectiveness:
                    best_effectiveness = effectiveness
                    best_choice_key = key
                    best_choice_new_area = new_area

        if best_choice_key is None:
            print("No further coverage improvement possible. Stopping.")
            break

        # Place the best scanner found
        scanner_type, px, py = best_choice_key
        cost = scanner_costs[scanner_type]
        
        total_cost += cost
        total_area += best_choice_new_area
        scanner_counts[scanner_type] += 1
        
        coverage_grid |= precomputed_masks[best_choice_key]

        print(f"Step {step}: Placed {scanner_type} at ({px}, {py}). "
              f"Coverage: {100*total_area/ROOM_AREA:.2f}%. Cost: ${total_cost}")
        step += 1

    # 5. Output Results
    print("\n--- Optimization Finished ---")
    print(f"Final covered area: {total_area} m^2")
    print(f"Final coverage: {100 * total_area / ROOM_AREA:.2f}% (Target was {COVERAGE_GOAL*100}%)")

    print("\nOptimal Scanner Configuration:")
    for scanner_type, count in scanner_counts.items():
        print(f"- {scanner_type}: {count} units")

    print("\nFinal Cost Calculation:")
    cost_str_parts = []
    for scanner_type, count in sorted(scanner_counts.items()):
        if count > 0:
            cost_str_parts.append(f"{count} * {scanner_costs[scanner_type]}")

    final_equation = f"{' + '.join(cost_str_parts)} = {total_cost}"
    print(final_equation)
    
    print(f"\n<<<The optimal total cost is {total_cost}>>>")


solve_scanner_placement()