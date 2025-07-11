import numpy as np
import math

def solve_scanner_placement():
    """
    Solves the museum scanner placement problem using a greedy algorithm.
    """
    # 1. Define Problem Constants
    WIDTH, HEIGHT = 140, 110
    TOTAL_AREA = WIDTH * HEIGHT
    TARGET_COVERAGE_RATIO = 0.88
    TARGET_COVERAGE_AREA = TOTAL_AREA * TARGET_COVERAGE_RATIO

    SCANNERS = {
        'C2': {'type': 'circle', 'radius': 20, 'cost': 20000},
        'C1': {'type': 'circle', 'radius': 5,  'cost': 1600},
        'R1': {'type': 'square', 'side': 10,   'cost': 2000}
    }

    print("--- Problem Setup ---")
    print(f"Room Area: {WIDTH}m x {HEIGHT}m = {TOTAL_AREA} m^2")
    print(f"Target Coverage: {TARGET_COVERAGE_RATIO:.0%} = {TARGET_COVERAGE_AREA:.2f} m^2")
    print("-" * 25)

    # 2. Setup Grids
    # Grid for placing scanner centers (coordinates are multiples of 5)
    placement_x_coords = range(0, WIDTH + 1, 5)
    placement_y_coords = range(0, HEIGHT + 1, 5)
    placement_points = [(x, y) for x in placement_x_coords for y in placement_y_coords]

    # Grid for tracking coverage (1 cell = 1 m^2). Size is +1 to include boundaries.
    coverage_grid = np.zeros((HEIGHT + 1, WIDTH + 1), dtype=bool)

    # Pre-compute coordinate grids for faster calculations
    yy, xx = np.mgrid[:HEIGHT + 1, :WIDTH + 1]

    # 3. Initialize State for Greedy Loop
    total_cost = 0
    total_covered_area = 0
    scanner_counts = {'C2': 0, 'C1': 0, 'R1': 0}
    
    print("Starting optimization... this may take a few minutes.\n")

    # 4. Main Greedy Loop
    while total_covered_area < TARGET_COVERAGE_AREA:
        best_value = -1
        best_placement = None
        best_mask = None

        # Iterate through all scanner types and all possible locations
        for name, props in SCANNERS.items():
            cost = props['cost']
            for cx, cy in placement_points:
                
                # Create a boolean mask representing the scanner's full coverage area
                if props['type'] == 'circle':
                    radius_sq = props['radius']**2
                    current_mask = (xx - cx)**2 + (yy - cy)**2 <= radius_sq
                else: # square
                    side = props['side']
                    half_side = side / 2
                    x_start = max(0, round(cx - half_side))
                    x_end = min(WIDTH + 1, round(cx + half_side))
                    y_start = max(0, round(cy - half_side))
                    y_end = min(HEIGHT + 1, round(cy + half_side))
                    current_mask = np.zeros_like(coverage_grid, dtype=bool)
                    current_mask[y_start:y_end, x_start:x_end] = True

                # Calculate new area (pixels on the mask that are NOT yet covered)
                newly_covered_mask = current_mask & ~coverage_grid
                new_area = np.sum(newly_covered_mask)

                if new_area == 0:
                    continue

                # Calculate value (cost-effectiveness for this specific placement)
                value = new_area / cost

                if value > best_value:
                    best_value = value
                    best_placement = (name, (cx, cy))
                    best_mask = newly_covered_mask
        
        if best_placement is None:
            print("Could not find any more effective placements. Stopping.")
            break

        # 5. Update State with the best choice found in this iteration
        scanner_name, location = best_placement
        
        # Update coverage grid by applying the new mask
        coverage_grid |= best_mask
        
        # Update totals
        total_cost += SCANNERS[scanner_name]['cost']
        newly_added_area = np.sum(best_mask)
        total_covered_area += newly_added_area
        scanner_counts[scanner_name] += 1

        print(f"Placed {scanner_name} @{location}. "
              f"Coverage: {total_covered_area / TOTAL_AREA:.2%}. Cost: ${total_cost}")

    # 6. Final Output
    print("\n--- Optimization Complete ---")

    c2_count = scanner_counts['C2']
    c1_count = scanner_counts['C1']
    r1_count = scanner_counts['R1']
    c2_cost = SCANNERS['C2']['cost']
    c1_cost = SCANNERS['C1']['cost']
    r1_cost = SCANNERS['R1']['cost']

    print(f"Final scanner configuration:")
    print(f"- C2 Scanners (20m radius): {c2_count}")
    print(f"- C1 Scanners (5m radius):  {c1_count}")
    print(f"- R1 Scanners (10m square): {r1_count}")
    print(f"\nFinal Coverage: {total_covered_area} / {TOTAL_AREA} = {total_covered_area / TOTAL_AREA:.4f} (Target was {TARGET_COVERAGE_RATIO:.4f})")
    
    print("\nOptimal Total Cost Calculation:")
    
    equation_parts = []
    if c2_count > 0:
        equation_parts.append(f"{c2_count} * {c2_cost}")
    if c1_count > 0:
        equation_parts.append(f"{c1_count} * {c1_cost}")
    if r1_count > 0:
        equation_parts.append(f"{r1_count} * {r1_cost}")

    final_equation = " + ".join(equation_parts)
    print(f"{final_equation} = {total_cost}")

    print(f"\n<<<{total_cost}>>>")

if __name__ == '__main__':
    solve_scanner_placement()