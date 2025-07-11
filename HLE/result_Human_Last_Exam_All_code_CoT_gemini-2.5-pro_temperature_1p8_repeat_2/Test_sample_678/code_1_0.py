import numpy as np
import math

def solve_scanner_placement():
    """
    Solves the scanner placement optimization problem using a greedy algorithm.
    """
    # --- Constants and Configuration ---
    ROOM_W, ROOM_H = 140, 110
    GRID_STEP = 5
    COVERAGE_TARGET_RATIO = 0.88

    TOTAL_AREA_POINTS = ROOM_W * ROOM_H
    TARGET_COVERAGE_POINTS = math.ceil(TOTAL_AREA_POINTS * COVERAGE_TARGET_RATIO)

    # Scanner properties: Cost, shape, and size (radius or half-side)
    SCANNERS_INFO = {
        'C2': {'cost': 20000, 'shape': 'circle', 'size': 20},
        'C1': {'cost': 1600, 'shape': 'circle', 'size': 5},
        'R1': {'cost': 2000, 'shape': 'square', 'size': 5}
    }

    # --- State Initialization ---
    # The coverage grid tracks which 1x1m squares are covered.
    # False = not covered, True = covered.
    coverage_grid = np.zeros((ROOM_H, ROOM_W), dtype=bool)

    # Create a set of all possible placement center points.
    placement_locs_x = np.arange(0, ROOM_W + 1, GRID_STEP)
    placement_locs_y = np.arange(0, ROOM_H + 1, GRID_STEP)
    available_placements = set()
    for x in placement_locs_x:
        for y in placement_locs_y:
            available_placements.add((x, y))

    total_cost = 0
    covered_points = 0
    placed_scanners_list = []

    # --- Pre-calculate scanner masks ---
    scanner_masks = {}
    for name, info in SCANNERS_INFO.items():
        size = info['size']
        dim = 2 * size + 1
        center = size
        yy, xx = np.mgrid[-center:center+1, -center:center+1]
        if info['shape'] == 'circle':
            mask = xx**2 + yy**2 <= size**2
        else: # square
            mask = np.ones((dim, dim), dtype=bool)
        scanner_masks[name] = mask

    # --- Main Greedy Loop ---
    print("Finding the optimal scanner placement. This may take a moment...")
    iteration = 0
    while covered_points < TARGET_COVERAGE_POINTS:
        iteration += 1
        best_choice = {'effectiveness': -1}

        # Iterate through each scanner type and each available location
        for name, info in SCANNERS_INFO.items():
            cost = info['cost']
            mask = scanner_masks[name]
            size = info['size']
            
            for loc in available_placements:
                center_x, center_y = loc

                # Define the slice of the room this scanner would affect
                y_min = max(0, center_y - size)
                y_max = min(ROOM_H, center_y + size + 1)
                x_min = max(0, center_x - size)
                x_max = min(ROOM_W, center_x + size + 1)
                
                # Clip the scanner mask to fit within the room boundaries
                mask_y_start = y_min - (center_y - size)
                mask_y_end = mask_y_start + (y_max - y_min)
                mask_x_start = x_min - (center_x - size)
                mask_x_end = mask_x_start + (x_max - x_min)
                clipped_mask = mask[mask_y_start:mask_y_end, mask_x_start:mask_x_end]

                # Get the corresponding slice from the main coverage grid
                room_slice = coverage_grid[y_min:y_max, x_min:x_max]
                
                # Calculate how many *new* points this scanner would cover
                newly_covered_mask = np.logical_and(clipped_mask, np.logical_not(room_slice))
                new_coverage_count = np.sum(newly_covered_mask)

                if new_coverage_count == 0:
                    continue

                effectiveness = new_coverage_count / cost

                if effectiveness > best_choice['effectiveness']:
                    best_choice.update({
                        'effectiveness': effectiveness,
                        'name': name,
                        'location': loc,
                        'new_coverage': new_coverage_count,
                        'cost': cost,
                        'placement_info': (y_min, y_max, x_min, x_max, clipped_mask)
                    })

        if best_choice['effectiveness'] == -1:
            print("Error: Cannot cover any more area. The target might be unreachable.")
            break
        
        # --- Apply the Best Choice Found in This Iteration ---
        bc = best_choice
        name, loc, cost = bc['name'], bc['location'], bc['cost']
        new_coverage, placement_info = bc['new_coverage'], bc['placement_info']

        # Update state
        total_cost += cost
        covered_points += new_coverage
        placed_scanners_list.append({'type': name, 'location': loc})
        available_placements.remove(loc)

        # Update the master coverage grid
        y_min, y_max, x_min, x_max, clipped_mask = placement_info
        room_slice_to_update = coverage_grid[y_min:y_max, x_min:x_max]
        np.logical_or(room_slice_to_update, clipped_mask, out=room_slice_to_update)
        
        print(f"  Step {iteration}: Placed {name} at {loc}. Coverage: {covered_points/TOTAL_AREA_POINTS:.2%}")

    # --- Final Results ---
    print("\n--- Optimization Complete ---")
    
    counts = {'C2': 0, 'C1': 0, 'R1': 0}
    for scanner in placed_scanners_list:
        counts[scanner['type']] += 1

    cost_c2 = SCANNERS_INFO['C2']['cost']
    cost_c1 = SCANNERS_INFO['C1']['cost']
    cost_r1 = SCANNERS_INFO['R1']['cost']

    print("Optimal Scanner Configuration:")
    print(f"- C2 Scanners (Cost: {cost_c2}): {counts['C2']}")
    print(f"- C1 Scanners (Cost: {cost_c1}): {counts['C1']}")
    print(f"- R1 Scanners (Cost: {cost_r1}): {counts['R1']}")
    print("\nPerformance:")
    print(f"- Final Coverage: {covered_points / TOTAL_AREA_POINTS:.4f} ({covered_points}/{TOTAL_AREA_POINTS} sq.m)")
    print(f"- Target Coverage: {COVERAGE_TARGET_RATIO:.4f} ({TARGET_COVERAGE_POINTS}/{TOTAL_AREA_POINTS} sq.m)")
    
    print("\nFinal Cost Calculation:")
    final_cost_eq = f"{counts['C2']} * {cost_c2} + {counts['C1']} * {cost_c1} + {counts['R1']} * {cost_r1} = {total_cost}"
    print(final_cost_eq)

    print(f"\nOptimal Total Cost: {total_cost}")
    print(f"\n<<<${total_cost}>>>")

if __name__ == '__main__':
    solve_scanner_placement()