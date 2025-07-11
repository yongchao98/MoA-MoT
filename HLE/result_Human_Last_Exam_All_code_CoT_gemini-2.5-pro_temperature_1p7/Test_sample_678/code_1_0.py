import numpy as np
import math

def solve_scanner_placement():
    """
    Solves the scanner placement optimization problem using a greedy algorithm.
    """
    # 1. Define constants and setup the environment
    ROOM_W = 110  # meters
    ROOM_L = 140  # meters
    ROOM_AREA = ROOM_W * ROOM_L
    TARGET_COVERAGE_RATIO = 0.88
    TARGET_COVERAGE_POINTS = math.floor(ROOM_AREA * TARGET_COVERAGE_RATIO)

    SCANNERS = {
        'C2': {'shape': 'circle', 'radius': 20, 'side': None, 'cost': 20000},
        'C1': {'shape': 'circle', 'radius': 5,  'side': None, 'cost': 1600}, # 10m diameter
        'R1': {'shape': 'square', 'side': 10, 'radius': None, 'cost': 2000}
    }

    # Grid for placing scanner centers
    X_CENTERS = range(0, ROOM_L + 1, 5)
    Y_CENTERS = range(0, ROOM_W + 1, 5)

    # High-resolution grid for coverage calculation (1x1m cells)
    y_coords, x_coords = np.mgrid[0:ROOM_W + 1, 0:ROOM_L + 1]
    TOTAL_GRID_POINTS = (ROOM_W + 1) * (ROOM_L + 1)

    # 2. Initialize variables for the greedy algorithm
    covered_mask = np.zeros((ROOM_W + 1, ROOM_L + 1), dtype=bool)
    total_cost = 0
    placements = []
    
    # Store counts of each scanner
    scanner_counts = {name: 0 for name in SCANNERS.keys()}

    print("Starting optimization... this may take a minute.")
    print(f"Target coverage: {TARGET_COVERAGE_POINTS} m^2 ({TARGET_COVERAGE_RATIO:.0%})")
    
    # 3. Main greedy loop
    iteration = 0
    while np.sum(covered_mask) < TARGET_COVERAGE_POINTS:
        iteration += 1
        best_option = {
            'metric': -1,
            'type': None,
            'center': None,
            'cost': float('inf'),
            'new_coverage': 0,
            'new_mask': None
        }

        # Iterate through every possible scanner placement to find the best one for this step
        for name, props in SCANNERS.items():
            for cx in X_CENTERS:
                for cy in Y_CENTERS:
                    # Generate the mask for this scanner
                    if props['shape'] == 'circle':
                        radius_sq = props['radius']**2
                        scanner_mask = (x_coords - cx)**2 + (y_coords - cy)**2 <= radius_sq
                    else:  # square
                        half_side = props['side'] / 2
                        scanner_mask = (np.abs(x_coords - cx) <= half_side) & (np.abs(y_coords - cy) <= half_side)

                    # Calculate newly covered points
                    newly_covered_mask = scanner_mask & ~covered_mask
                    new_coverage_points = np.sum(newly_covered_mask)
                    
                    if new_coverage_points == 0:
                        continue
                        
                    metric = new_coverage_points / props['cost']
                    
                    if metric > best_option['metric']:
                        best_option = {
                            'metric': metric,
                            'type': name,
                            'center': (cx, cy),
                            'cost': props['cost'],
                            'new_coverage': new_coverage_points,
                            'new_mask': newly_covered_mask
                        }

        # If no option adds coverage, stop to prevent an infinite loop
        if best_option['metric'] == -1:
            print("Could not find any more placements to add coverage. Stopping.")
            break

        # Add the best scanner found in this iteration
        total_cost += best_option['cost']
        covered_mask |= best_option['new_mask']
        scanner_counts[best_option['type']] += 1

        current_coverage_area = np.sum(covered_mask)
        coverage_ratio = current_coverage_area / TOTAL_GRID_POINTS
        
        print(f"Step {iteration}: Added {best_option['type']} at ({best_option['center'][0]}m, {best_option['center'][1]}m). Cost: {total_cost}. Coverage: {coverage_ratio:.2%}")

    # 4. Final report
    print("\n--- Optimization Complete ---")
    final_coverage_area = np.sum(covered_mask)
    final_coverage_ratio = final_coverage_area / TOTAL_GRID_POINTS
    
    print(f"Target Coverage: {TARGET_COVERAGE_POINTS} m^2 ({TARGET_COVERAGE_RATIO:.2%})")
    print(f"Achieved Coverage: {final_coverage_area} m^2 ({final_coverage_ratio:.2%})")
    
    print("\nFinal Scanner Configuration:")
    print(f"- Type C2 (20m radius, cost 20000): {scanner_counts['C2']}")
    print(f"- Type C1 (5m radius, cost 1600):  {scanner_counts['C1']}")
    print(f"- Type R1 (10m square, cost 2000): {scanner_counts['R1']}")
    
    # Format the cost calculation equation
    cost_components = []
    if scanner_counts['C2'] > 0:
        cost_components.append(f"{scanner_counts['C2']} * {SCANNERS['C2']['cost']}")
    if scanner_counts['C1'] > 0:
        cost_components.append(f"{scanner_counts['C1']} * {SCANNERS['C1']['cost']}")
    if scanner_counts['R1'] > 0:
        cost_components.append(f"{scanner_counts['R1']} * {SCANNERS['R1']['cost']}")

    print("\nFinal Cost Calculation:")
    cost_equation = " + ".join(cost_components)
    print(f"{cost_equation} = {total_cost}")
    
    print(f"\nOptimal Total Cost: {total_cost}")
    return total_cost

if __name__ == '__main__':
    final_cost = solve_scanner_placement()
    print(f"\n<<< {final_cost} >>>")