import numpy as np
import math

def solve_scanner_placement():
    """
    Solves the museum scanner placement problem using a greedy algorithm.
    """
    # --- Step 1: Problem Definition & Setup ---
    ROOM_W = 140  # meters
    ROOM_H = 110  # meters
    TOTAL_AREA = ROOM_W * ROOM_H
    COVERAGE_TARGET_RATIO = 0.88
    COVERAGE_TARGET_AREA = TOTAL_AREA * COVERAGE_TARGET_RATIO
    PLACEMENT_STEP = 5

    scanners_def = [
        {'name': 'C2', 'shape': 'circle', 'radius': 20, 'cost': 20000},
        {'name': 'C1', 'shape': 'circle', 'radius': 5, 'cost': 1600},
        {'name': 'R1', 'shape': 'square', 'side': 10, 'cost': 2000}
    ]

    print("--- Problem Setup ---")
    print(f"Room Area: {TOTAL_AREA} sq.m.")
    print(f"Target Coverage: {COVERAGE_TARGET_RATIO*100}% ({COVERAGE_TARGET_AREA:.0f} sq.m.)")
    print("-" * 25)

    # --- Step 2: Pre-calculate all possible scanner coverage masks ---
    # Create a 1x1m grid representing the room floor
    yy, xx = np.mgrid[0:ROOM_H, 0:ROOM_W]

    # Generate all possible placement locations (centers are multiples of 5)
    possible_locs_x = range(0, ROOM_W + 1, PLACEMENT_STEP)
    possible_locs_y = range(0, ROOM_H + 1, PLACEMENT_STEP)
    placement_locations = [(x, y) for x in possible_locs_x for y in possible_locs_y]

    print("Pre-calculating all possible scanner coverage masks...")
    candidate_scanners = []
    for loc in placement_locations:
        cx, cy = loc
        for scanner_type in scanners_def:
            if scanner_type['shape'] == 'circle':
                r = scanner_type['radius']
                mask = (xx - cx)**2 + (yy - cy)**2 <= r**2
            elif scanner_type['shape'] == 'square':
                s = scanner_type['side']
                mask = (np.abs(xx - cx) <= s / 2) & (np.abs(yy - cy) <= s / 2)
            
            # Only add candidates that cover some area within the room
            if np.any(mask):
                candidate_scanners.append({
                    'name': scanner_type['name'],
                    'location': loc,
                    'cost': scanner_type['cost'],
                    'mask': mask
                })
    print(f"Generated {len(candidate_scanners)} candidate scanner placements.")
    print("-" * 25)


    # --- Step 3: Run the greedy algorithm ---
    print("Starting greedy optimization...\n")
    coverage_grid = np.full((ROOM_H, ROOM_W), False, dtype=bool)
    total_cost = 0
    solution_counts = {s['name']: 0 for s in scanners_def}
    iteration = 0

    while np.sum(coverage_grid) < COVERAGE_TARGET_AREA:
        iteration += 1
        best_candidate = None
        max_efficiency = -1
        
        for candidate in candidate_scanners:
            # Calculate the new area this candidate would cover
            newly_covered_mask = candidate['mask'] & ~coverage_grid
            newly_covered_area = np.sum(newly_covered_mask)
            
            if newly_covered_area > 0:
                efficiency = newly_covered_area / candidate['cost']
                if efficiency > max_efficiency:
                    max_efficiency = efficiency
                    best_candidate = candidate

        if best_candidate is None:
            print("No more coverage can be added. Stopping.")
            break
            
        # Add the best found scanner to the solution
        total_cost += best_candidate['cost']
        solution_counts[best_candidate['name']] += 1
        coverage_grid |= best_candidate['mask']
        
        current_coverage_area = np.sum(coverage_grid)
        current_coverage_ratio = current_coverage_area / TOTAL_AREA

        print(f"Iter {iteration}: Added {best_candidate['name']} at {best_candidate['location']}. "
              f"Coverage: {current_coverage_ratio:.4f}. Total Cost: {total_cost}")

        # Remove all candidates at the chosen location
        chosen_location = best_candidate['location']
        candidate_scanners = [c for c in candidate_scanners if c['location'] != chosen_location]

    # --- Step 4: Report the final results ---
    print("\n" + "-" * 25)
    print("--- Optimization Finished ---")
    final_coverage_area = np.sum(coverage_grid)
    final_coverage_ratio = final_coverage_area / TOTAL_AREA

    print(f"\nFinal Coverage Achieved: {final_coverage_ratio:.4f} ({final_coverage_area} sq.m.)")
    print("\nOptimal Scanner Configuration:")
    for name, count in sorted(solution_counts.items()):
        if count > 0:
            print(f"- {name}: {count} units")

    print("\nFinal Equation:")
    cost_str_parts = []
    for name, count in sorted(solution_counts.items()):
        if count > 0:
            cost = next(s['cost'] for s in scanners_def if s['name'] == name)
            cost_str_parts.append(f"{count} * {cost}")
    
    # Using print to ensure the equation and final answer are clearly displayed
    final_equation_str = f"Total Cost = {' + '.join(cost_str_parts)} = {total_cost}"
    print(final_equation_str)
    
    # Final answer in the specified format
    print(f"\n<<<The optimal total cost is {total_cost}>>>")


if __name__ == '__main__':
    solve_scanner_placement()