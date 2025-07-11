import numpy as np
import math

def solve_scanner_placement():
    """
    Solves the museum scanner placement optimization problem using a greedy algorithm.
    """
    # --- Problem Definition ---

    # Room dimensions and area
    ROOM_WIDTH = 140
    ROOM_HEIGHT = 110
    ROOM_AREA = ROOM_WIDTH * ROOM_HEIGHT

    # Coverage requirement
    COVERAGE_RATIO = 0.88
    # The target area is based on discrete 1x1m grid points
    TARGET_COVERAGE_AREA = ROOM_AREA * COVERAGE_RATIO

    # Scanner definitions
    SCANNERS = [
        {
            'name': 'C2',
            'type': 'circle',
            'param': 20,  # radius
            'cost': 20000
        },
        {
            'name': 'C1',
            'type': 'circle',
            'param': 5,   # radius (from 10m diameter)
            'cost': 1600
        },
        {
            'name': 'R1',
            'type': 'square',
            'param': 10,  # side length
            'cost': 2000
        }
    ]

    # Possible placement coordinates (multiples of 5m)
    STEP = 5
    placement_coords_x = np.arange(0, ROOM_WIDTH + 1, STEP)
    placement_coords_y = np.arange(0, ROOM_HEIGHT + 1, STEP)

    # --- Optimization using a Greedy Algorithm ---

    # Discretize the room area into a grid of 1x1m points.
    # We use +1 because dimensions are inclusive (e.g., 0m to 140m).
    grid_x, grid_y = np.mgrid[0:ROOM_WIDTH + 1, 0:ROOM_HEIGHT + 1]
    
    # State-tracking variables
    coverage_grid = np.zeros((ROOM_WIDTH + 1, ROOM_HEIGHT + 1), dtype=bool)
    total_cost = 0
    placed_scanners_count = {'C2': 0, 'C1': 0, 'R1': 0}

    # Loop until the coverage target is met
    while np.sum(coverage_grid) < TARGET_COVERAGE_AREA:
        best_effectiveness = -1
        best_candidate = None

        # Iterate through all scanner types and all possible placements
        # to find the most cost-effective one to add next.
        for scanner in SCANNERS:
            for cx in placement_coords_x:
                for cy in placement_coords_y:
                    # Generate a mask representing the scanner's coverage area
                    if scanner['type'] == 'circle':
                        radius = scanner['param']
                        mask = (grid_x - cx)**2 + (grid_y - cy)**2 <= radius**2
                    elif scanner['type'] == 'square':
                        side = scanner['param']
                        half_side = side / 2
                        mask = (np.abs(grid_x - cx) <= half_side) & (np.abs(grid_y - cy) <= half_side)
                    
                    # Calculate area that is newly covered (respecting overlaps)
                    newly_covered_mask = mask & ~coverage_grid
                    newly_covered_area = np.sum(newly_covered_mask)

                    if newly_covered_area == 0:
                        continue

                    # The "bang for the buck"
                    effectiveness = newly_covered_area / scanner['cost']

                    if effectiveness > best_effectiveness:
                        best_effectiveness = effectiveness
                        best_candidate = {
                            'name': scanner['name'],
                            'cost': scanner['cost'],
                            'mask': newly_covered_mask,
                        }

        if best_candidate is None:
            # This case happens if no scanner can add any new coverage
            print("Could not add more coverage to meet the target. Halting.")
            break

        # Add the best scanner found in this iteration
        total_cost += best_candidate['cost']
        coverage_grid |= best_candidate['mask']
        placed_scanners_count[best_candidate['name']] += 1

    # --- Final Output ---
    print("Optimization Result:")

    num_c2 = placed_scanners_count['C2']
    num_c1 = placed_scanners_count['C1']
    num_r1 = placed_scanners_count['R1']
    cost_c2 = SCANNERS[0]['cost']
    cost_c1 = SCANNERS[1]['cost']
    cost_r1 = SCANNERS[2]['cost']

    print("\nThe optimal cost is calculated by the following equation:")
    # We print each number as requested
    print(f"{num_c2} * {cost_c2} + {num_c1} * {cost_c1} + {num_r1} * {cost_r1} = {total_cost}")

    final_coverage_pts = np.sum(coverage_grid)
    # Total points in the grid representation
    total_possible_pts = (ROOM_WIDTH + 1) * (ROOM_HEIGHT + 1)
    final_coverage_ratio = final_coverage_pts / total_possible_pts

    print(f"\nThis configuration achieves {final_coverage_ratio:.2%} coverage, meeting the {COVERAGE_RATIO:.0%} target.")
    print(f"The optimal total cost found is: {total_cost}")

    print(f"\n<<<{total_cost}>>>")


if __name__ == '__main__':
    solve_scanner_placement()