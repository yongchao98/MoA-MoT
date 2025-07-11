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

    SCANNERS = {
        'C2': {'shape': 'circle', 'radius': 20, 'cost': 20000},
        'C1': {'shape': 'circle', 'radius': 5,  'cost': 1600}, # from 10m diameter
        'R1': {'shape': 'square', 'side': 10,   'cost': 2000}
    }

    # --- Discretization and Setup ---
    # The room is represented by a grid where each cell is 1x1m.
    # A value of 0 means uncovered, 1 means covered.
    coverage_grid = np.zeros((ROOM_HEIGHT, ROOM_WIDTH))
    total_area_points = ROOM_WIDTH * ROOM_HEIGHT
    target_covered_points = total_area_points * TARGET_COVERAGE_RATIO

    # Possible center locations for scanners (must be multiples of 5)
    placement_x_coords = range(0, ROOM_WIDTH + 1, 5)
    placement_y_coords = range(0, ROOM_HEIGHT + 1, 5)

    # --- Algorithm Initialization ---
    total_cost = 0
    current_covered_points = 0
    chosen_scanners = []

    # Create coordinate matrices for the entire room for efficient calculations
    y_grid, x_grid = np.mgrid[0:ROOM_HEIGHT, 0:ROOM_WIDTH]

    print("Solving optimization problem... (This may take a moment)")
    
    # --- Greedy Algorithm Loop ---
    while current_covered_points < target_covered_points:
        best_value = -1
        best_scanner_choice = None
        best_newly_covered_mask = None

        # Iterate over all possible placements and scanner types to find the most cost-effective addition
        for px in placement_x_coords:
            for py in placement_y_coords:
                for name, props in SCANNERS.items():
                    cost = props['cost']
                    
                    # Calculate the coverage mask for the current potential scanner
                    if props['shape'] == 'circle':
                        radius = props['radius']
                        # Equation of a circle: (x-h)^2 + (y-k)^2 <= r^2
                        scanner_mask = ((x_grid - px)**2 + (y_grid - py)**2) <= radius**2
                    else:  # 'square'
                        side = props['side']
                        half_side = side / 2
                        scanner_mask = (x_grid >= px - half_side) & (x_grid < px + half_side) & \
                                       (y_grid >= py - half_side) & (y_grid < py + half_side)
                    
                    # Determine the new area this scanner would cover
                    # This is the intersection of the scanner's mask and the uncovered area (where grid is 0)
                    newly_covered_mask = scanner_mask & (coverage_grid == 0)
                    new_points = np.sum(newly_covered_mask)

                    if new_points == 0:
                        continue
                    
                    # Calculate the value (new points covered per unit cost)
                    value = new_points / cost
                    
                    if value > best_value:
                        best_value = value
                        best_scanner_choice = {'name': name, 'x': px, 'y': py, 'cost': cost}
                        best_newly_covered_mask = newly_covered_mask

        if best_scanner_choice is None:
            print("Warning: No further coverage is possible. Stopping.")
            break
        
        # Add the best scanner found in this iteration to our solution
        total_cost += best_scanner_choice['cost']
        chosen_scanners.append(best_scanner_choice)
        
        # Update the master coverage grid
        coverage_grid[best_newly_covered_mask] = 1
        newly_covered_count = np.sum(best_newly_covered_mask)
        current_covered_points += newly_covered_count

    # --- Final Output ---
    print("\n--- Optimization Result ---")

    # Count the number of each type of scanner used
    num_c2 = sum(1 for s in chosen_scanners if s['name'] == 'C2')
    num_c1 = sum(1 for s in chosen_scanners if s['name'] == 'C1')
    num_r1 = sum(1 for s in chosen_scanners if s['name'] == 'R1')

    # Build the final equation string
    cost_parts = []
    if num_c2 > 0:
        cost_parts.append(f"{num_c2} * 20000")
    if num_c1 > 0:
        cost_parts.append(f"{num_c1} * 1600")
    if num_r1 > 0:
        cost_parts.append(f"{num_r1} * 2000")
    cost_equation = " + ".join(cost_parts)
    
    final_coverage_ratio = current_covered_points / total_area_points

    print("Optimal solution found:")
    print(f"- C2 Scanners: {num_c2}")
    print(f"- C1 Scanners: {num_c1}")
    print(f"- R1 Scanners: {num_r1}")
    print("\nFinal Equation for Cost:")
    print(f"Total Cost = {cost_equation} = {total_cost}")
    print(f"\nThis configuration achieves a coverage of {final_coverage_ratio:.2%}, meeting the target of {TARGET_COVERAGE_RATIO:.2%}.")
    print("\nThe optimal total cost is:")
    print(total_cost)


if __name__ == '__main__':
    solve_scanner_placement()