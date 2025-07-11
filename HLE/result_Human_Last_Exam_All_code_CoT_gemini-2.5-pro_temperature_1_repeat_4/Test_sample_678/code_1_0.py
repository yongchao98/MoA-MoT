import math

def solve_scanner_placement():
    """
    Solves the museum scanner placement problem using a greedy algorithm.
    """
    # 1. Define problem constants
    ROOM_WIDTH = 140  # meters
    ROOM_HEIGHT = 110 # meters
    GRID_STEP = 5     # meters for scanner placement
    COVERAGE_GOAL = 0.88

    SCANNERS_SPECS = {
        'C2': {'shape': 'circle', 'radius': 20, 'cost': 20000},
        'C1': {'shape': 'circle', 'radius': 5,  'cost': 1600}, # 10m diameter
        'R1': {'shape': 'square', 'side': 10,   'cost': 2000}
    }

    TOTAL_CELLS = ROOM_WIDTH * ROOM_HEIGHT
    TARGET_COVERAGE_CELLS = math.floor(TOTAL_CELLS * COVERAGE_GOAL)

    # 2. Pre-compute all potential scanner installations
    def get_covered_cells(spec, center_x, center_y):
        """Calculates the set of 1x1m cells covered by a scanner."""
        cells = set()
        if spec['shape'] == 'circle':
            r = spec['radius']
            # Iterate through the bounding box of the circle for efficiency
            for x in range(max(0, math.floor(center_x - r)), min(ROOM_WIDTH, math.ceil(center_x + r))):
                for y in range(max(0, math.floor(center_y - r)), min(ROOM_HEIGHT, math.ceil(center_y + r))):
                    # Check if the center of the cell is within the circle
                    if (x + 0.5 - center_x)**2 + (y + 0.5 - center_y)**2 <= r**2:
                        cells.add((x, y))
        elif spec['shape'] == 'square':
            half_side = spec['side'] / 2
            min_x = max(0, int(center_x - half_side))
            max_x = min(ROOM_WIDTH, int(center_x + half_side))
            min_y = max(0, int(center_y - half_side))
            max_y = min(ROOM_HEIGHT, int(center_y + half_side))
            for x in range(min_x, max_x):
                for y in range(min_y, max_y):
                    cells.add((x, y))
        return cells

    potential_scanners = []
    possible_locations = []
    for x in range(0, ROOM_WIDTH + 1, GRID_STEP):
        for y in range(0, ROOM_HEIGHT + 1, GRID_STEP):
            possible_locations.append((x, y))

    for name, spec in SCANNERS_SPECS.items():
        for loc in possible_locations:
            cells = get_covered_cells(spec, loc[0], loc[1])
            if cells: # Only add if it covers any area
                potential_scanners.append({
                    'name': name,
                    'location': loc,
                    'cost': spec['cost'],
                    'cells': cells
                })

    # 3. Execute the greedy algorithm
    total_cost = 0
    covered_cells_set = set()
    solution_counts = {'C2': 0, 'C1': 0, 'R1': 0}
    
    # Keep track of used locations to enforce one scanner per location
    used_locations = set()

    while len(covered_cells_set) < TARGET_COVERAGE_CELLS:
        best_scanner = None
        max_efficiency = -1
        best_new_cells = set()

        # Find the most cost-efficient scanner to add
        for scanner in potential_scanners:
            if scanner['location'] in used_locations:
                continue

            newly_covered = scanner['cells'] - covered_cells_set
            
            if not newly_covered:
                continue
            
            efficiency = len(newly_covered) / scanner['cost']
            
            if efficiency > max_efficiency:
                max_efficiency = efficiency
                best_scanner = scanner
                best_new_cells = newly_covered

        if best_scanner is None:
            # This happens if no remaining scanner can add new coverage
            print("Stopping: No more coverage can be added.")
            break

        # Add the best scanner found to the solution
        name = best_scanner['name']
        cost = best_scanner['cost']
        location = best_scanner['location']

        total_cost += cost
        covered_cells_set.update(best_new_cells)
        solution_counts[name] += 1
        used_locations.add(location)

    # 4. Print the final results
    cost_c2 = SCANNERS_SPECS['C2']['cost']
    cost_c1 = SCANNERS_SPECS['C1']['cost']
    cost_r1 = SCANNERS_SPECS['R1']['cost']
    
    num_c2 = solution_counts['C2']
    num_c1 = solution_counts['C1']
    num_r1 = solution_counts['R1']

    print("Optimal Scanner Configuration and Cost:")
    print(f"{num_c2} * {cost_c2} + {num_c1} * {cost_c1} + {num_r1} * {cost_r1} = {total_cost}")
    
    actual_coverage_percent = len(covered_cells_set) / TOTAL_CELLS * 100
    print(f"\nFinal Coverage: {len(covered_cells_set)} / {TOTAL_CELLS} cells = {actual_coverage_percent:.2f}%")
    print(f"Target Coverage: {TARGET_COVERAGE_CELLS} cells ({COVERAGE_GOAL * 100}%)")
    print(f"Final Total Cost: {total_cost}")
    
    print(f"<<<{total_cost}>>>")


if __name__ == '__main__':
    solve_scanner_placement()
