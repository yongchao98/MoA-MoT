import numpy as np
import math

def solve_scanner_placement():
    """
    Solves the scanner placement optimization problem using a greedy algorithm.
    """
    # --- Problem Configuration ---
    ROOM_W = 140  # Room width in meters
    ROOM_H = 110  # Room height in meters
    TARGET_COVERAGE_RATIO = 0.88
    CENTER_SPACING = 5 # Scanners centers are multiples of 5m

    SCANNERS = [
        {'name': 'C2', 'cost': 20000, 'type': 'circle', 'param': 20}, # radius 20m
        {'name': 'C1', 'cost': 1600, 'type': 'circle', 'param': 5},   # radius 5m (diameter 10m)
        {'name': 'R1', 'cost': 2000, 'type': 'square', 'param': 10},  # side 10m
    ]

    # --- Algorithm Setup ---
    TOTAL_AREA = ROOM_W * ROOM_H
    TARGET_AREA = TOTAL_AREA * TARGET_COVERAGE_RATIO

    print(f"Room Dimensions: {ROOM_W}m x {ROOM_H}m (Total Area: {TOTAL_AREA} m^2)")
    print(f"Target Coverage: {TARGET_COVERAGE_RATIO:.0%} ({TARGET_AREA} m^2)\n")

    # Generate a grid of all possible scanner center locations
    possible_centers_x = range(0, ROOM_W + 1, CENTER_SPACING)
    possible_centers_y = range(0, ROOM_H + 1, CENTER_SPACING)
    possible_centers = [(x, y) for x in possible_centers_x for y in possible_centers_y]

    # Create a high-resolution grid representing the room for coverage calculation.
    # Each cell is 1m x 1m. We check the center of each cell for coverage.
    yy, xx = np.mgrid[0:ROOM_H, 0:ROOM_W]
    pixel_centers_x = xx + 0.5
    pixel_centers_y = yy + 0.5

    # --- State Variables ---
    coverage_grid = np.zeros((ROOM_H, ROOM_W), dtype=bool)
    total_cost = 0
    scanner_counts = {s['name']: 0 for s in SCANNERS}
    
    # --- Main Greedy Loop ---
    iteration = 0
    while np.sum(coverage_grid) < TARGET_AREA:
        iteration += 1
        best_option = {
            'marginal_cost': float('inf'),
            'scanner': None,
            'center': None,
            'mask': None,
        }

        # Find the most cost-effective scanner to add in this step
        for scanner in SCANNERS:
            for center in possible_centers:
                cx, cy = center
                
                # Generate a boolean mask for the area this scanner would cover
                if scanner['type'] == 'circle':
                    radius_sq = scanner['param']**2
                    mask = (pixel_centers_x - cx)**2 + (pixel_centers_y - cy)**2 <= radius_sq
                elif scanner['type'] == 'square':
                    half_side = scanner['param'] / 2
                    mask = (np.abs(pixel_centers_x - cx) <= half_side) & \
                           (np.abs(pixel_centers_y - cy) <= half_side)
                
                # Calculate the newly covered area by this scanner
                newly_covered_area = np.sum(mask & ~coverage_grid)

                if newly_covered_area > 0:
                    marginal_cost = scanner['cost'] / newly_covered_area
                    if marginal_cost < best_option['marginal_cost']:
                        best_option['marginal_cost'] = marginal_cost
                        best_option['scanner'] = scanner
                        best_option['center'] = center
                        best_option['mask'] = mask

        # If no more area can be covered, break the loop
        if best_option['scanner'] is None:
            print("Could not find any more area to cover. Stopping.")
            break

        # Add the best scanner found in this iteration
        chosen_scanner = best_option['scanner']
        total_cost += chosen_scanner['cost']
        scanner_counts[chosen_scanner['name']] += 1
        coverage_grid |= best_option['mask'] # Update the master coverage grid
        
        current_covered_area = np.sum(coverage_grid)
        coverage_percent = (current_covered_area / TOTAL_AREA) * 100
        
        print(f"Step {iteration}: Added {chosen_scanner['name']} at {best_option['center']}. "
              f"Total Cost: {total_cost}, Coverage: {coverage_percent:.2f}%")

    # --- Final Output ---
    print("\n--- Optimization Complete ---")
    
    final_covered_area = np.sum(coverage_grid)
    final_coverage_ratio = final_covered_area / TOTAL_AREA
    print(f"Final Achieved Coverage: {final_covered_area} m^2 ({final_coverage_ratio:.2%})")
    
    print("\nFinal Scanner Configuration:")
    for name, count in scanner_counts.items():
        print(f"- {name}: {count} units")

    print("\nFinal Cost Equation:")
    cost_details = []
    for scanner in SCANNERS:
        name = scanner['name']
        cost = scanner['cost']
        count = scanner_counts[name]
        if count > 0:
            cost_details.append(f"{count} * {cost}")
    
    # Ensure all scanner types are in the equation, even if count is 0
    if scanner_counts['C2'] == 0: cost_details.insert(0, f"0 * {SCANNERS[0]['cost']}")
    if scanner_counts['C1'] == 0: cost_details.insert(1, f"0 * {SCANNERS[1]['cost']}")
    if scanner_counts['R1'] == 0: cost_details.insert(2, f"0 * {SCANNERS[2]['cost']}")


    print(f"{' + '.join(cost_details)} = {total_cost}")
    
    print(f"\nOptimal Total Cost: {total_cost}")
    print(f"<<<{total_cost}>>>")

if __name__ == '__main__':
    solve_scanner_placement()