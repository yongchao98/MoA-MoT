import math
import numpy as np

def solve_scanner_placement():
    """
    Solves the museum scanner placement problem using a greedy algorithm.
    """
    # --- Problem Definition ---
    ROOM_WIDTH = 140
    ROOM_HEIGHT = 110
    TARGET_COVERAGE_RATIO = 0.88
    
    SCANNERS = {
        'C2': {'shape': 'circle', 'radius': 20, 'cost': 20000},
        'C1': {'shape': 'circle', 'radius': 5, 'cost': 1600},  # diameter 10m
        'R1': {'shape': 'square', 'side': 10, 'cost': 2000},
    }
    
    # Grid for placing scanner centers
    CENTER_STEP = 5
    # Grid for approximating area coverage
    COVERAGE_STEP = 1

    # 1. Model the Room with a grid of points
    points_x = np.arange(0, ROOM_WIDTH + COVERAGE_STEP, COVERAGE_STEP)
    points_y = np.arange(0, ROOM_HEIGHT + COVERAGE_STEP, COVERAGE_STEP)
    total_points_in_room = len(points_x) * len(points_y)
    target_covered_points_count = math.ceil(total_points_in_room * TARGET_COVERAGE_RATIO)

    # 2. Define all possible scanner placements
    center_coords_x = np.arange(0, ROOM_WIDTH + CENTER_STEP, CENTER_STEP)
    center_coords_y = np.arange(0, ROOM_HEIGHT + CENTER_STEP, CENTER_STEP)
    
    available_placements = {} # Using dict for easier removal of chosen locations
    for cx in center_coords_x:
        for cy in center_coords_y:
            center_key = (cx, cy)
            available_placements[center_key] = []
            for name, props in SCANNERS.items():
                placement = {
                    'name': name,
                    'center': (cx, cy),
                    'cost': props['cost'],
                    'props': props,
                    'key': (name, (cx, cy))
                }
                available_placements[center_key].append(placement)

    # 3. Pre-calculate coverage for each placement
    coverage_map = {}
    for center_key, placements_at_center in available_placements.items():
        for p in placements_at_center:
            covered_points = set()
            cx, cy = p['center']
            props = p['props']
            
            if props['shape'] == 'circle':
                r = props['radius']
                r_sq = r**2
                # Bounding box for efficiency
                min_x, max_x = cx - r, cx + r
                min_y, max_y = cy - r, cy + r
            elif props['shape'] == 'square':
                s_half = props['side'] / 2
                min_x, max_x = cx - s_half, cx + s_half
                min_y, max_y = cy - s_half, cy + s_half
            
            # Filter points within the scanner's bounding box to speed up check
            points_to_check_x = points_x[(points_x >= min_x) & (points_x <= max_x)]
            points_to_check_y = points_y[(points_y >= min_y) & (points_y <= max_y)]
            
            for x in points_to_check_x:
                for y in points_to_check_y:
                    is_covered = False
                    if props['shape'] == 'circle':
                        if (x - cx)**2 + (y - cy)**2 <= r_sq:
                            is_covered = True
                    elif props['shape'] == 'square':
                        if abs(x - cx) <= props['side'] / 2 and abs(y - cy) <= props['side'] / 2:
                            is_covered = True
                    if is_covered:
                        covered_points.add((x, y))
            coverage_map[p['key']] = covered_points

    # 4. Greedy Selection Algorithm
    total_cost = 0
    total_covered_points_set = set()
    chosen_scanners = []

    while len(total_covered_points_set) < target_covered_points_count:
        best_effectiveness = -1
        best_placement = None
        
        # Find the most cost-effective scanner among all remaining placements
        locations_to_check = list(available_placements.keys())
        for center_key in locations_to_check:
            for p in available_placements[center_key]:
                points_for_this_scanner = coverage_map[p['key']]
                newly_covered_points_count = len(points_for_this_scanner - total_covered_points_set)
                
                if newly_covered_points_count == 0:
                    continue
                
                effectiveness = newly_covered_points_count / p['cost']
                if effectiveness > best_effectiveness:
                    best_effectiveness = effectiveness
                    best_placement = p

        if best_placement is None:
            # No scanner can add new coverage, so we stop.
            break

        # "Place" the best scanner found
        total_cost += best_placement['cost']
        total_covered_points_set.update(coverage_map[best_placement['key']])
        chosen_scanners.append(best_placement)
        
        # Remove the chosen location from the pool of available placements
        del available_placements[best_placement['center']]

    # 5. Output the result
    counts = {'C2': 0, 'C1': 0, 'R1': 0}
    for scanner in chosen_scanners:
        counts[scanner['name']] += 1

    final_coverage_ratio = len(total_covered_points_set) / total_points_in_room

    print("--- Optimization Result ---")
    print(f"Target Coverage: {TARGET_COVERAGE_RATIO:.2%} | Achieved Coverage: {final_coverage_ratio:.2%}\n")
    print("Optimal Scanner Configuration:")
    print(f"- C2 Scanners (cost 20000 each): {counts['C2']}")
    print(f"- C1 Scanners (cost 1600 each): {counts['C1']}")
    print(f"- R1 Scanners (cost 2000 each): {counts['R1']}\n")
    print("Final cost calculation:")
    
    # Build and print the equation string
    cost_c2 = counts['C2'] * SCANNERS['C2']['cost']
    cost_c1 = counts['C1'] * SCANNERS['C1']['cost']
    cost_r1 = counts['R1'] * SCANNERS['R1']['cost']
    
    print(f"{counts['C2']} * {SCANNERS['C2']['cost']} + {counts['C1']} * {SCANNERS['C1']['cost']} + {counts['R1']} * {SCANNERS['R1']['cost']} = {total_cost}")

    print(f"\n<<< {total_cost} >>>")


solve_scanner_placement()