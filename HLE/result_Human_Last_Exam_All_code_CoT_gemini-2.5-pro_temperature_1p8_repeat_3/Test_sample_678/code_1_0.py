import math
from collections import Counter

def solve_scanner_placement():
    """
    Solves the museum scanner placement problem by finding a cost-effective
    configuration that meets the coverage target.
    """
    # --- Step 1: Define Problem Constants ---
    ROOM_WIDTH = 140
    ROOM_HEIGHT = 110
    TARGET_COVERAGE_RATIO = 0.88
    PLACEMENT_STEP = 5  # Scanner centers must be multiples of 5m

    SCANNERS = {
        'C2': {'shape': 'circle', 'radius': 20, 'cost': 20000},
        'C1': {'shape': 'circle', 'radius': 5, 'cost': 1600}, # diameter 10m
        'R1': {'shape': 'square', 'side': 10, 'cost': 2000},
    }

    # --- Step 2: Discretize the Room into a Grid ---
    # A smaller grid step increases accuracy but is much slower. 2.5m is a good balance.
    GRID_STEP = 2.5
    
    grid_points = set()
    for i in range(int(ROOM_WIDTH / GRID_STEP)):
        for j in range(int(ROOM_HEIGHT / GRID_STEP)):
            grid_points.add((i * GRID_STEP, j * GRID_STEP))
            
    total_points = len(grid_points)
    target_points_to_cover = int(total_points * TARGET_COVERAGE_RATIO)
    
    # --- Step 3: Pre-compute Coverage for All Possible Scanners ---
    possible_locations = []
    for x in range(0, ROOM_WIDTH + 1, PLACEMENT_STEP):
        for y in range(0, ROOM_HEIGHT + 1, PLACEMENT_STEP):
            possible_locations.append((x, y))

    potential_scanners = {}
    for name, spec in SCANNERS.items():
        for loc in possible_locations:
            cx, cy = loc
            covered = set()
            # Determine the bounding box to check only relevant points
            if spec['shape'] == 'circle':
                radius = spec['radius']
                radius_sq = radius ** 2
                for px, py in grid_points:
                    if (px - cx)**2 + (py - cy)**2 < radius_sq:
                        covered.add((px, py))
            elif spec['shape'] == 'square':
                half_side = spec['side'] / 2
                for px, py in grid_points:
                    if (cx - half_side <= px < cx + half_side) and \
                       (cy - half_side <= py < cy + half_side):
                        covered.add((px, py))

            if covered:
                potential_scanners[(name, loc)] = {
                    'cost': spec['cost'],
                    'coverage': covered
                }

    # --- Step 4: Run Greedy Algorithm to Select Scanners ---
    covered_points = set()
    selected_scanners = []
    total_cost = 0

    while len(covered_points) < target_points_to_cover:
        best_scanner_key = None
        best_effectiveness = 0  # Maximize (new_points / cost)
        best_new_points = set()

        for key, scanner_data in potential_scanners.items():
            newly_covered = scanner_data['coverage'].difference(covered_points)
            if not newly_covered:
                continue

            effectiveness = len(newly_covered) / scanner_data['cost']
            if effectiveness > best_effectiveness:
                best_effectiveness = effectiveness
                best_scanner_key = key
                best_new_points = newly_covered

        if best_scanner_key is None:
            # This case happens if no scanner can add new coverage
            break
            
        name, loc = best_scanner_key
        cost = SCANNERS[name]['cost']
        
        selected_scanners.append({'type': name, 'cost': cost})
        total_cost += cost
        covered_points.update(best_new_points)

    # --- Step 5: Print the Final Results ---
    final_coverage = len(covered_points) / total_points
    scanner_counts = Counter(s['type'] for s in selected_scanners)

    print("--- Optimization Result ---")
    print(f"Optimal Total Cost: {total_cost}")
    print(f"Final Coverage Achieved: {final_coverage:.2%}\n")
    
    cost_equation_parts = []
    print("Cost breakdown:")
    for scanner_type, count in sorted(scanner_counts.items()):
        cost = SCANNERS[scanner_type]['cost']
        cost_equation_parts.append(f"{count} * {cost}")
    
    final_equation = " + ".join(cost_equation_parts)
    print(f"Final equation: {final_equation} = {total_cost}")

    # Print each number in the final equation as requested
    print("\nEach number in the final equation:")
    for part in cost_equation_parts:
        # Split by ' * ' and print numbers
        nums = part.split(' * ')
        print(nums[0])
        print(nums[1])

    # Final answer format
    print(f"<<<{total_cost}>>>")

if __name__ == '__main__':
    solve_scanner_placement()