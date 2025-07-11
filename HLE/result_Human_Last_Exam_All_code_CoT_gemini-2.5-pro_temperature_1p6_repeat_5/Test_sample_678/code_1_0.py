import numpy as np
import math

def solve_scanner_placement():
    """
    Formulates the scanner placement problem and solves it using a greedy algorithm
    to find a cost-effective solution for covering a designated area.
    """
    # --- Configuration ---
    ROOM_WIDTH = 140  # meters
    ROOM_HEIGHT = 110 # meters
    COVERAGE_TARGET_FRACTION = 0.88
    STEP = 5 # meters, for placement grid

    SCANNERS = [
        {
            "name": "C2",
            "type": "circle",
            "param": 20,  # radius
            "cost": 20000
        },
        {
            "name": "C1",
            "type": "circle",
            "param": 5,   # radius (from 10m diameter)
            "cost": 1600
        },
        {
            "name": "R1",
            "type": "square",
            "param": 10,  # side
            "cost": 2000
        }
    ]

    # --- Setup ---
    TOTAL_AREA_POINTS = ROOM_WIDTH * ROOM_HEIGHT
    TARGET_COVERAGE_POINTS = TOTAL_AREA_POINTS * COVERAGE_TARGET_FRACTION

    # Discretize the room into a 1x1m grid.
    # We use (height, width) convention for numpy arrays, which corresponds to (y, x)
    coverage_grid = np.zeros((ROOM_HEIGHT, ROOM_WIDTH), dtype=bool)

    # Generate possible placement locations (centers are multiples of 5m)
    placement_x = range(0, ROOM_WIDTH + 1, STEP)
    placement_y = range(0, ROOM_HEIGHT + 1, STEP)
    available_placements = [(x, y) for x in placement_x for y in placement_y]

    total_cost = 0
    placed_scanners_count = {"C2": 0, "C1": 0, "R1": 0}

    # --- Main Greedy Loop ---
    current_coverage_points = 0
    while current_coverage_points < TARGET_COVERAGE_POINTS:
        best_choice = None
        max_value = -1.0  # cost-effectiveness (area per cost)

        # Iterate through every possible scanner at every available location
        for scanner_info in SCANNERS:
            cost = scanner_info["cost"]
            for center_x, center_y in available_placements:

                # Determine scanner bounding box, clamped to room dimensions
                if scanner_info["type"] == "circle":
                    radius = scanner_info["param"]
                    y_min = max(0, center_y - radius)
                    y_max = min(ROOM_HEIGHT, center_y + radius)
                    x_min = max(0, center_x - radius)
                    x_max = min(ROOM_WIDTH, center_x + radius)
                else:  # square
                    half_side = scanner_info["param"] // 2
                    y_min = max(0, center_y - half_side)
                    y_max = min(ROOM_HEIGHT, center_y + half_side)
                    x_min = max(0, center_x - half_side)
                    x_max = min(ROOM_WIDTH, center_x + half_side)

                if y_min >= y_max or x_min >= x_max:
                    continue

                # Get the sub-grid slice of uncovered points
                sub_grid_slice = (slice(y_min, y_max), slice(x_min, x_max))
                uncovered_sub_grid = ~coverage_grid[sub_grid_slice]

                # Create a mask for the scanner shape within the bounding box
                yy, xx = np.ogrid[y_min:y_max, x_min:x_max]
                if scanner_info["type"] == "circle":
                    radius = scanner_info["param"]
                    mask = (xx + 0.5 - center_x)**2 + (yy + 0.5 - center_y)**2 < radius**2
                else: # square
                    half_side = scanner_info["param"] / 2.0
                    mask = (np.abs(xx + 0.5 - center_x) < half_side) & (np.abs(yy + 0.5 - center_y) < half_side)

                # Calculate new area by applying the mask to the uncovered points
                newly_covered_area = np.sum(mask & uncovered_sub_grid)

                if newly_covered_area > 0:
                    value = newly_covered_area / cost
                    if value > max_value:
                        max_value = value
                        best_choice = {
                            "scanner": scanner_info,
                            "location": (center_x, center_y),
                            "grid_slice": sub_grid_slice,
                            "mask_to_apply": mask,
                        }

        if best_choice is None:
            print("No further coverage possible with remaining placements. Stopping.")
            break
        
        # --- Place the Best Scanner Found ---
        scanner = best_choice["scanner"]
        location = best_choice["location"]
        
        total_cost += scanner["cost"]
        placed_scanners_count[scanner["name"]] += 1
        
        # Update coverage grid efficiently
        grid_slice = best_choice["grid_slice"]
        mask = best_choice["mask_to_apply"]
        coverage_grid[grid_slice][mask] = True

        current_coverage_points = np.sum(coverage_grid)
        available_placements.remove(location)

    # --- Final Output ---
    print("Optimization Complete.")
    print(f"Target coverage: {TARGET_COVERAGE_POINTS} m^2 ({COVERAGE_TARGET_FRACTION:.2%})")
    print(f"Achieved coverage: {current_coverage_points} m^2 ({current_coverage_points / TOTAL_AREA_POINTS:.2%})")
    print("\n--- Optimal Scanner Configuration ---")
    
    cost_parts = []
    for name, count in placed_scanners_count.items():
        if count > 0:
            cost = next(s['cost'] for s in SCANNERS if s['name'] == name)
            print(f"Number of {name} scanners: {count}")
            cost_parts.append(f"{count} * {cost}")
    
    print("\nFinal Equation for Total Cost:")
    final_equation = " + ".join(cost_parts) + f" = {total_cost}"
    print(final_equation)
    print(f"<<<{int(total_cost)}>>>")

if __name__ == '__main__':
    solve_scanner_placement()