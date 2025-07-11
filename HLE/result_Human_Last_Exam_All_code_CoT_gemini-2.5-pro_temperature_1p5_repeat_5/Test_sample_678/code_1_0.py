import numpy as np
import math

def solve_scanner_placement():
    """
    Solves the museum scanner placement problem using a greedy algorithm.
    """
    # --- Problem Definition ---
    ROOM_WIDTH = 140
    ROOM_HEIGHT = 110
    TARGET_COVERAGE_RATIO = 0.88
    
    ROOM_AREA = ROOM_WIDTH * ROOM_HEIGHT
    TARGET_COVERED_AREA = ROOM_AREA * TARGET_COVERAGE_RATIO

    scanners_def = {
        "C2": {"shape": "circle", "radius": 20, "cost": 20000},
        "C1": {"shape": "circle", "radius": 5, "cost": 1600},
        "R1": {"shape": "square", "side": 10, "cost": 2000}
    }

    # --- Grid Setup ---
    # A 1x1m grid to represent the room area for coverage calculation
    coverage_grid = np.zeros((ROOM_WIDTH, ROOM_HEIGHT), dtype=bool)

    # Possible center locations for scanners (multiples of 5m)
    placement_step = 5
    possible_x_coords = range(0, ROOM_WIDTH + 1, placement_step)
    possible_y_coords = range(0, ROOM_HEIGHT + 1, placement_step)
    
    # --- Greedy Algorithm ---
    total_cost = 0
    covered_area = 0
    selected_scanners_log = []

    print("Running optimization... (this may take a moment)")
    
    iteration = 0
    while covered_area < TARGET_COVERED_AREA:
        iteration += 1
        best_value = -1
        best_scanner_choice = None

        # Iterate through all scanner types and all possible locations
        for name, scanner_props in scanners_def.items():
            for cx in possible_x_coords:
                for cy in possible_y_coords:
                    
                    # Calculate how many *new* cells this scanner would cover
                    newly_covered_count = 0
                    if scanner_props['shape'] == 'circle':
                        radius = scanner_props['radius']
                        r_squared = radius**2
                        min_x = max(0, math.floor(cx - radius))
                        max_x = min(ROOM_WIDTH, math.ceil(cx + radius))
                        min_y = max(0, math.floor(cy - radius))
                        max_y = min(ROOM_HEIGHT, math.ceil(cy + radius))
                        for x in range(min_x, max_x):
                            for y in range(min_y, max_y):
                                if not coverage_grid[x, y]:
                                    if (x + 0.5 - cx)**2 + (y + 0.5 - cy)**2 <= r_squared:
                                        newly_covered_count += 1
                    
                    elif scanner_props['shape'] == 'square':
                        half_side = scanner_props['side'] / 2
                        min_x = max(0, math.floor(cx - half_side))
                        max_x = min(ROOM_WIDTH, math.ceil(cx + half_side))
                        min_y = max(0, math.floor(cy - half_side))
                        max_y = min(ROOM_HEIGHT, math.ceil(cy + half_side))
                        for x in range(min_x, max_x):
                            for y in range(min_y, max_y):
                                if not coverage_grid[x, y]:
                                    newly_covered_count += 1
                    
                    if newly_covered_count > 0:
                        cost = scanner_props['cost']
                        # Value = New Coverage per unit cost
                        value = newly_covered_count / cost
                        
                        if value > best_value:
                            best_value = value
                            best_scanner_choice = {
                                "name": name,
                                "center": (cx, cy),
                                "cost": cost,
                                "new_coverage": newly_covered_count
                            }
        
        # If no placement can add new coverage, stop.
        if best_scanner_choice is None:
            print("Could not add more coverage. Stopping.")
            break
        
        # Add the best found scanner to our solution
        choice_name = best_scanner_choice["name"]
        choice_center = best_scanner_choice["center"]
        choice_cost = best_scanner_choice["cost"]
        choice_new_coverage = best_scanner_choice["new_coverage"]

        selected_scanners_log.append(best_scanner_choice)
        total_cost += choice_cost
        covered_area += choice_new_coverage

        # Update the coverage grid based on the chosen scanner
        props = scanners_def[choice_name]
        if props['shape'] == 'circle':
            radius = props['radius']
            r_squared = radius**2
            min_x = max(0, math.floor(choice_center[0] - radius))
            max_x = min(ROOM_WIDTH, math.ceil(choice_center[0] + radius))
            min_y = max(0, math.floor(choice_center[1] - radius))
            max_y = min(ROOM_HEIGHT, math.ceil(choice_center[1] + radius))
            for x in range(min_x, max_x):
                for y in range(min_y, max_y):
                    if (x + 0.5 - choice_center[0])**2 + (y + 0.5 - choice_center[1])**2 <= r_squared:
                        coverage_grid[x, y] = True
        elif props['shape'] == 'square':
            half_side = props['side'] / 2
            min_x = max(0, math.floor(choice_center[0] - half_side))
            max_x = min(ROOM_WIDTH, math.ceil(choice_center[0] + half_side))
            min_y = max(0, math.floor(choice_center[1] - half_side))
            max_y = min(ROOM_HEIGHT, math.ceil(choice_center[1] + half_side))
            for x in range(min_x, max_x):
                for y in range(min_y, max_y):
                    coverage_grid[x, y] = True

        print(f"Step {iteration}: Added {choice_name} @ {choice_center}. "
              f"Coverage: {covered_area / ROOM_AREA:.2%}. Cost: {total_cost}")

    # --- Final Analysis & Output ---
    print("\n--- Optimization Complete ---")
    final_coverage_ratio = covered_area / ROOM_AREA
    print(f"Target Coverage: {TARGET_COVERAGE_RATIO:.2%} ({int(TARGET_COVERED_AREA)} sq m)")
    print(f"Achieved Coverage: {final_coverage_ratio:.2%} ({covered_area} sq m)")

    # Count the number of each scanner type
    scanner_counts = {name: 0 for name in scanners_def.keys()}
    for scanner in selected_scanners_log:
        scanner_counts[scanner['name']] += 1

    print("\n--- Final Cost Calculation ---")
    
    cost_components = []
    for name, count in scanner_counts.items():
        if count > 0:
            cost = scanners_def[name]['cost']
            cost_components.append(f"{count} * {cost}")
    
    equation = " + ".join(cost_components)
    print(f"Optimal scanner configuration cost:")
    print(f"{equation} = {total_cost}")
    
    print("\n--- Final Answer ---")
    print(f"<<<{total_cost}>>>")

if __name__ == '__main__':
    solve_scanner_placement()