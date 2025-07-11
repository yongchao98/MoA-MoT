import numpy as np
import math

def solve_scanner_placement():
    """
    Formulates and solves the scanner placement optimization problem using a greedy algorithm.
    """
    # 1. Define problem parameters
    ROOM_WIDTH = 140
    ROOM_HEIGHT = 110
    TARGET_COVERAGE_RATIO = 0.88
    PLACEMENT_GRID_STEP = 5

    # Scanner types: (name, cost, shape, parameter)
    SCANNER_TYPES = [
        {'name': 'C2', 'cost': 20000, 'shape': 'circle', 'param': 20}, # radius
        {'name': 'C1', 'cost': 1600,  'shape': 'circle', 'param': 5},  # radius (from 10m diameter)
        {'name': 'R1', 'cost': 2000,  'shape': 'square', 'param': 10}, # side length
    ]

    # 2. Set up the space and grids
    # Grid of points to check for coverage (1m resolution)
    y_coords, x_coords = np.mgrid[0:ROOM_HEIGHT + 1, 0:ROOM_WIDTH + 1]
    total_points = (ROOM_WIDTH + 1) * (ROOM_HEIGHT + 1)
    target_covered_points = math.ceil(total_points * TARGET_COVERAGE_RATIO)

    # Potential scanner center locations (multiples of 5)
    x_centers = list(range(0, ROOM_WIDTH + 1, PLACEMENT_GRID_STEP))
    y_centers = list(range(0, ROOM_HEIGHT + 1, PLACEMENT_GRID_STEP))
    possible_locations = [(x, y) for x in x_centers for y in y_centers]

    # 3. Initialize the greedy algorithm
    coverage_grid = np.zeros_like(x_coords, dtype=bool)
    current_covered_points = 0
    total_cost = 0
    placed_scanners = []

    print("Starting optimization... please wait.")

    # 4. Run the greedy algorithm loop
    while current_covered_points < target_covered_points:
        best_effectiveness = -1
        best_scanner_info = None
        best_scanner_mask = None
        
        # Iterate through all possible scanner placements to find the most cost-effective one
        for location in possible_locations:
            for scanner_type in SCANNER_TYPES:
                cx, cy = location
                
                # Calculate the coverage mask for this scanner candidate
                if scanner_type['shape'] == 'circle':
                    r = scanner_type['param']
                    mask = (x_coords - cx)**2 + (y_coords - cy)**2 <= r**2
                elif scanner_type['shape'] == 'square':
                    s_half = scanner_type['param'] / 2
                    mask = ((x_coords >= cx - s_half) & (x_coords <= cx + s_half) &
                            (y_coords >= cy - s_half) & (y_coords <= cy + s_half))
                
                # Calculate how many *new* points this scanner would cover
                newly_covered_mask = np.logical_and(mask, np.logical_not(coverage_grid))
                new_points_count = np.sum(newly_covered_mask)

                if new_points_count == 0:
                    continue

                cost = scanner_type['cost']
                effectiveness = new_points_count / cost

                if effectiveness > best_effectiveness:
                    best_effectiveness = effectiveness
                    best_scanner_info = {
                        'type': scanner_type['name'],
                        'center': location,
                        'cost': cost
                    }
                    best_scanner_mask = newly_covered_mask

        # If no effective scanner can be found, stop.
        if best_scanner_info is None:
            print("Warning: Could not find any more effective scanners to place.")
            break

        # Add the best scanner found in this iteration to the solution
        placed_scanners.append(best_scanner_info)
        total_cost += best_scanner_info['cost']
        coverage_grid = np.logical_or(coverage_grid, best_scanner_mask)
        current_covered_points += np.sum(best_scanner_mask)
        
    # 5. Print the final results
    print("\nOptimization Complete.")
    print("="*30)
    
    final_coverage = np.sum(coverage_grid)
    print(f"Target Coverage: >= {TARGET_COVERAGE_RATIO * 100}% ({target_covered_points} points)")
    print(f"Achieved Coverage: {final_coverage / total_points * 100:.2f}% ({final_coverage} points)")
    
    num_c2 = sum(1 for s in placed_scanners if s['type'] == 'C2')
    num_c1 = sum(1 for s in placed_scanners if s['type'] == 'C1')
    num_r1 = sum(1 for s in placed_scanners if s['type'] == 'R1')
    
    print("\nOptimal Scanner Configuration:")
    print(f"- C2 Scanners (20m radius, cost 20000): {num_c2}")
    print(f"- C1 Scanners (10m diameter, cost 1600): {num_c1}")
    print(f"- R1 Scanners (10m side, cost 2000): {num_r1}")
    
    print("\nFinal Cost Calculation:")
    cost_c2 = 20000
    cost_c1 = 1600
    cost_r1 = 2000
    equation = f"{num_c2} * {cost_c2} + {num_c1} * {cost_c1} + {num_r1} * {cost_r1} = {total_cost}"
    print(equation)
    
    print("\nOptimal Total Cost:")
    print(total_cost)
    
    return total_cost

if __name__ == '__main__':
    final_cost = solve_scanner_placement()
    print(f"\n<<<ANSWERS>>>")
    print(f"<<<{final_cost}>>>")
