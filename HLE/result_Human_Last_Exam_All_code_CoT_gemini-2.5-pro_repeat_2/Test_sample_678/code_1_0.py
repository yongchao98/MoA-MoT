import numpy as np
import math

def solve_scanner_placement():
    """
    Solves the museum scanner placement problem using a greedy algorithm.
    """
    # --- Configuration ---
    ROOM_WIDTH = 140
    ROOM_HEIGHT = 110
    GRID_SPACING = 5
    TARGET_COVERAGE_RATIO = 0.88

    SCANNERS = {
        'C2': {'shape': 'circle', 'param': 20, 'cost': 20000},  # param is radius
        'C1': {'shape': 'circle', 'param': 5,  'cost': 1600},   # param is radius (from 10m diameter)
        'R1': {'shape': 'square', 'param': 10, 'cost': 2000}    # param is side length
    }

    # --- Initialization ---
    total_room_area = ROOM_WIDTH * ROOM_HEIGHT
    target_area = total_room_area * TARGET_COVERAGE_RATIO

    # Create a grid of points representing the center of each 1m x 1m square
    yy, xx = np.mgrid[0:ROOM_HEIGHT, 0:ROOM_WIDTH]
    point_centers_y = yy + 0.5
    point_centers_x = xx + 0.5

    # Grid for placing scanners (centers must be multiples of 5)
    placement_xs = np.arange(0, ROOM_WIDTH + 1, GRID_SPACING)
    placement_ys = np.arange(0, ROOM_HEIGHT + 1, GRID_SPACING)
    placement_locs = [(x, y) for x in placement_xs for y in placement_ys]

    # State variables for the greedy algorithm
    covered_grid = np.zeros((ROOM_HEIGHT, ROOM_WIDTH), dtype=bool)
    total_cost = 0
    total_coverage = 0
    placed_scanners = []

    # --- Greedy Algorithm Loop ---
    while total_coverage < target_area:
        best_candidate = {
            'ratio': -1,
            'name': None,
            'center': None,
            'cost': -1,
            'new_coverage_mask': None,
            'new_area': 0
        }

        # Find the most cost-effective scanner to add in this step
        for name, specs in SCANNERS.items():
            cost = specs['cost']
            if cost <= 0: continue

            for center_x, center_y in placement_locs:
                # Calculate the coverage mask for this scanner candidate
                if specs['shape'] == 'circle':
                    radius = specs['param']
                    # Check if the center of each grid cell is within the circle
                    mask = (point_centers_x - center_x)**2 + (point_centers_y - center_y)**2 <= radius**2
                elif specs['shape'] == 'square':
                    half_side = specs['param'] / 2
                    # Check if the center of each grid cell is within the square
                    mask = (np.abs(point_centers_x - center_x) <= half_side) & \
                           (np.abs(point_centers_y - center_y) <= half_side)

                # Calculate new area this scanner would cover by ignoring already covered cells
                newly_covered_mask = mask & (~covered_grid)
                new_area = np.sum(newly_covered_mask)

                if new_area > 0:
                    ratio = new_area / cost
                    if ratio > best_candidate['ratio']:
                        best_candidate.update({
                            'ratio': ratio,
                            'name': name,
                            'center': (center_x, center_y),
                            'cost': cost,
                            'new_coverage_mask': newly_covered_mask,
                            'new_area': new_area
                        })

        # If no placement can add new coverage, stop.
        if best_candidate['name'] is None:
            print("Warning: Could not add more coverage. Target may not be reached.")
            break

        # Add the best found scanner to our solution and update the state
        covered_grid |= best_candidate['new_coverage_mask']
        total_cost += best_candidate['cost']
        total_coverage += best_candidate['new_area']
        placed_scanners.append(best_candidate)

    # --- Final Results ---
    scanner_counts = {name: 0 for name in SCANNERS}
    for scanner in placed_scanners:
        scanner_counts[scanner['name']] += 1

    num_c2 = scanner_counts.get('C2', 0)
    num_c1 = scanner_counts.get('C1', 0)
    num_r1 = scanner_counts.get('R1', 0)

    print("Optimal solution found:")
    print(f"Achieved coverage: {total_coverage / total_room_area:.2%} ({int(total_coverage)} / {int(target_area)} sqm required)")
    print(f"Number of C2 scanners: {num_c2}")
    print(f"Number of C1 scanners: {num_c1}")
    print(f"Number of R1 scanners: {num_r1}")
    print("\nFinal cost equation:")
    print(f"20000 * {num_c2} + 1600 * {num_c1} + 2000 * {num_r1} = {total_cost}")
    print(f"\n<<<${total_cost}>>>")

if __name__ == '__main__':
    solve_scanner_placement()