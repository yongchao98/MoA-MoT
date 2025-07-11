import math
import numpy as np
import sys

def solve_scanner_placement():
    """
    Solves the museum scanner placement problem using a greedy algorithm.
    """
    # 1. Define Constants and Problem Space
    ROOM_WIDTH = 140
    ROOM_HEIGHT = 110
    PLACEMENT_STEP = 5
    TARGET_COVERAGE_RATIO = 0.88

    # Using a 1x1m grid for coverage calculation
    coverage_grid = np.zeros((ROOM_WIDTH + 1, ROOM_HEIGHT + 1), dtype=bool)
    total_points = (ROOM_WIDTH + 1) * (ROOM_HEIGHT + 1)
    target_covered_points = math.ceil(total_points * TARGET_COVERAGE_RATIO)

    # 2. Define Scanners
    scanners = [
        {'name': 'C2', 'cost': 20000, 'type': 'circle', 'radius': 20},
        {'name': 'C1', 'cost': 1600, 'type': 'circle', 'radius': 5},  # diameter 10m -> radius 5m
        {'name': 'R1', 'cost': 2000, 'type': 'square', 'side': 10}
    ]

    # 3. Define Placement Locations
    placement_locations = []
    for x in range(0, ROOM_WIDTH + 1, PLACEMENT_STEP):
        for y in range(0, ROOM_HEIGHT + 1, PLACEMENT_STEP):
            placement_locations.append((x, y))
            
    # Helper function to get newly covered points for a potential scanner
    def get_newly_covered_points_coords(center_x, center_y, scanner_info, grid):
        new_coords = []
        if scanner_info['type'] == 'circle':
            radius = scanner_info['radius']
            r_sq = radius ** 2
            min_px = max(0, center_x - radius)
            max_px = min(ROOM_WIDTH, center_x + radius)
            min_py = max(0, center_y - radius)
            max_py = min(ROOM_HEIGHT, center_y + radius)
            
            for px in range(min_px, max_px + 1):
                for py in range(min_py, max_py + 1):
                    if not grid[px, py]:
                        if (px - center_x)**2 + (py - center_y)**2 <= r_sq:
                            new_coords.append((px, py))
                            
        elif scanner_info['type'] == 'square':
            side = scanner_info['side']
            half_side = side / 2
            min_px = max(0, math.ceil(center_x - half_side))
            max_px = min(ROOM_WIDTH, math.floor(center_x + half_side))
            min_py = max(0, math.ceil(center_y - half_side))
            max_py = min(ROOM_HEIGHT, math.floor(center_y + half_side))
            
            for px in range(min_px, max_px + 1):
                for py in range(min_py, max_py + 1):
                    if not grid[px, py]:
                        new_coords.append((px, py))
        return new_coords

    # 4. Main Greedy Algorithm
    total_cost = 0
    solution_scanners = []
    current_covered_points = 0
    
    print("Starting optimization... this may take a moment.")
    sys.stdout.flush()

    while current_covered_points < target_covered_points:
        best_option = None
        max_effectiveness = -1.0

        for loc_x, loc_y in placement_locations:
            for scanner in scanners:
                newly_covered_coords = get_newly_covered_points_coords(loc_x, loc_y, scanner, coverage_grid)
                newly_covered_count = len(newly_covered_coords)

                if newly_covered_count > 0:
                    effectiveness = newly_covered_count / scanner['cost']
                    if effectiveness > max_effectiveness:
                        max_effectiveness = effectiveness
                        best_option = {
                            'scanner': scanner,
                            'new_points_coords': newly_covered_coords
                        }
        
        if best_option is None:
            print("Warning: Cannot find any placement to increase coverage. Stopping.")
            break

        # Add the best found scanner to our solution
        solution_scanners.append(best_option['scanner'])
        total_cost += best_option['scanner']['cost']
        new_points_count = len(best_option['new_points_coords'])
        current_covered_points += new_points_count
        
        # Update coverage grid
        for px, py in best_option['new_points_coords']:
            coverage_grid[px, py] = True
            
        progress_percent = 100 * current_covered_points / target_covered_points
        print(f"Progress: {progress_percent:.1f}% | Coverage: {100 * current_covered_points / total_points:.2f}% | Cost: {total_cost}", end='\r')


    # 5. Output Results
    print("\n\n--- Optimization Complete ---")

    num_c2 = sum(1 for s in solution_scanners if s['name'] == 'C2')
    num_c1 = sum(1 for s in solution_scanners if s['name'] == 'C1')
    num_r1 = sum(1 for s in solution_scanners if s['name'] == 'R1')

    final_coverage_ratio = current_covered_points / total_points
    cost_c2 = scanners[0]['cost']
    cost_c1 = scanners[1]['cost']
    cost_r1 = scanners[2]['cost']

    print(f"Optimal scanner configuration found:")
    print(f"- Type C2 (20m radius, {cost_c2} cost): {num_c2} units")
    print(f"- Type C1 (5m radius, {cost_c1} cost):  {num_c1} units")
    print(f"- Type R1 (10m side, {cost_r1} cost):   {num_r1} units")
    print(f"\nFinal coverage: {final_coverage_ratio:.4f} (Target was >= {TARGET_COVERAGE_RATIO})")
    print(f"Total points covered: {current_covered_points} / {total_points}")
    print("\nFinal cost calculation:")
    print(f"{num_c2} * {cost_c2} + {num_c1} * {cost_c1} + {num_r1} * {cost_r1} = {total_cost}")

    print(f"\nOptimal total cost: {total_cost}")
    return total_cost

if __name__ == '__main__':
    final_cost = solve_scanner_placement()
    print(f"\n<<<${final_cost}>>>")
