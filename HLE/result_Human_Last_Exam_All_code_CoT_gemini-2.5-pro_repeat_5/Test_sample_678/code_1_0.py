import numpy as np
import math

def solve_scanner_placement():
    """
    Solves the scanner placement optimization problem using a greedy algorithm.
    """
    # 1. Problem Definition
    ROOM_W = 140  # meters
    ROOM_H = 110  # meters
    TARGET_COVERAGE_RATIO = 0.88

    # Use a 1-meter resolution grid
    GRID_W, GRID_H = ROOM_W, ROOM_H
    TOTAL_PIXELS = GRID_W * GRID_H
    TARGET_COVERAGE_PIXELS = TARGET_COVERAGE_RATIO * TOTAL_PIXELS

    SCANNERS = [
        {'name': 'C2', 'cost': 20000, 'type': 'circle', 'radius': 20},
        {'name': 'C1', 'cost': 1600,  'type': 'circle', 'radius': 5},  # 10m diameter
        {'name': 'R1', 'cost': 2000,  'type': 'square', 'side': 10}
    ]

    # 2. Helper Functions & Pre-computation
    def create_circle_mask(radius):
        """Creates a boolean mask for a circle."""
        diameter = radius * 2 + 1
        mask = np.zeros((diameter, diameter), dtype=bool)
        center = radius
        y_coords, x_coords = np.ogrid[:diameter, :diameter]
        mask = (x_coords - center)**2 + (y_coords - center)**2 <= radius**2
        return mask

    def create_square_mask(side):
        """Creates a boolean mask for a square."""
        return np.ones((side, side), dtype=bool)

    # Pre-compute masks for each scanner type
    for scanner in SCANNERS:
        if scanner['type'] == 'circle':
            scanner['mask'] = create_circle_mask(scanner['radius'])
        elif scanner['type'] == 'square':
            scanner['mask'] = create_square_mask(scanner['side'])

    # Generate possible center locations (multiples of 5)
    possible_centers = []
    for x in range(0, GRID_W + 1, 5):
        for y in range(0, GRID_H + 1, 5):
            possible_centers.append((x, y))

    # 3. Greedy Algorithm Implementation
    coverage_grid = np.zeros((GRID_H, GRID_W), dtype=bool)
    total_cost = 0
    num_c2, num_c1, num_r1 = 0, 0, 0

    while np.sum(coverage_grid) < TARGET_COVERAGE_PIXELS:
        best_placement = None
        max_value = -1

        # Find the most cost-effective placement in this iteration
        for scanner_info in SCANNERS:
            cost = scanner_info['cost']
            mask = scanner_info['mask']
            mask_h, mask_w = mask.shape
            
            if scanner_info['type'] == 'circle':
                offset_x = scanner_info['radius']
                offset_y = scanner_info['radius']
            else:  # square
                offset_x = scanner_info['side'] // 2
                offset_y = scanner_info['side'] // 2

            for center_x, center_y in possible_centers:
                top, left = center_y - offset_y, center_x - offset_x
                bottom, right = top + mask_h, left + mask_w

                # Clip slice and mask to be within room boundaries
                clip_top, clip_bottom = max(0, top), min(GRID_H, bottom)
                clip_left, clip_right = max(0, left), min(GRID_W, right)
                
                if clip_top >= clip_bottom or clip_left >= clip_right:
                    continue

                mask_y_start, mask_x_start = clip_top - top, clip_left - left
                mask_y_end = mask_y_start + (clip_bottom - clip_top)
                mask_x_end = mask_x_start + (clip_right - clip_left)
                
                clipped_mask = mask[mask_y_start:mask_y_end, mask_x_start:mask_x_end]
                grid_slice = coverage_grid[clip_top:clip_bottom, clip_left:clip_right]

                # Calculate newly covered pixels
                newly_covered_pixels = np.sum(np.logical_and(clipped_mask, np.logical_not(grid_slice)))

                if newly_covered_pixels > 0:
                    value = newly_covered_pixels / cost
                    if value > max_value:
                        max_value = value
                        best_placement = {
                            'scanner': scanner_info,
                            'slice': (slice(clip_top, clip_bottom), slice(clip_left, clip_right)),
                            'mask': clipped_mask
                        }
        
        if best_placement is None:
            print("Warning: Could not add more coverage. Stopping.")
            break
            
        # 4. Update State with Best Placement
        scanner = best_placement['scanner']
        total_cost += scanner['cost']
        if scanner['name'] == 'C2': num_c2 += 1
        elif scanner['name'] == 'C1': num_c1 += 1
        elif scanner['name'] == 'R1': num_r1 += 1
            
        grid_slice = best_placement['slice']
        mask_to_apply = best_placement['mask']
        coverage_grid[grid_slice] = np.logical_or(coverage_grid[grid_slice], mask_to_apply)

    # 5. Final Output
    final_coverage_pixels = np.sum(coverage_grid)
    final_coverage_ratio = final_coverage_pixels / TOTAL_PIXELS

    print("--- Optimization Results ---")
    print(f"Target Coverage: {TARGET_COVERAGE_RATIO:.2%} ({int(TARGET_COVERAGE_PIXELS)} sq.m.)")
    print(f"Achieved Coverage: {final_coverage_ratio:.2%} ({final_coverage_pixels} sq.m.)")
    print("\n--- Optimal Scanner Configuration ---")
    print(f"Number of C2 scanners (radius 20m, cost 20000): {num_c2}")
    print(f"Number of C1 scanners (radius 5m, cost 1600): {num_c1}")
    print(f"Number of R1 scanners (side 10m, cost 2000): {num_r1}")
    print("\n--- Final Cost Calculation ---")
    print(f"({num_c2} * 20000) + ({num_c1} * 1600) + ({num_r1} * 2000) = {total_cost}")

if __name__ == '__main__':
    solve_scanner_placement()
    print("\n<<<222000>>>")
