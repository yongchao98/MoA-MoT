import math

def solve_scanner_optimization():
    """
    This script solves the scanner placement problem to achieve at least 88% coverage
    of a 140x110m room with the minimum cost using a greedy algorithm.
    """
    # --- Problem Definition ---
    ROOM_WIDTH = 140
    ROOM_HEIGHT = 110
    COVERAGE_TARGET = 0.88
    GRID_SPACING = 5  # For scanner centers

    SCANNER_TYPES = {
        'C2': {'shape': 'circle', 'radius': 20, 'cost': 20000},
        'C1': {'shape': 'circle', 'radius': 5, 'cost': 1600}, # 10m diameter
        'R1': {'shape': 'square', 'side': 10, 'cost': 2000}
    }

    # 1. Create grid of checkpoints for coverage calculation (1x1m grid)
    print("Step 1: Modeling the room with a 1x1m checkpoint grid...")
    checkpoints = []
    for x in range(ROOM_WIDTH + 1):
        for y in range(ROOM_HEIGHT + 1):
            checkpoints.append((x, y))
    total_points = len(checkpoints)
    target_covered_points = math.ceil(total_points * COVERAGE_TARGET)
    print(f"Room area is represented by {total_points} points. Target is to cover {target_covered_points} points.")

    # 2. Identify all possible scanner locations
    scanner_locations = []
    for x in range(0, ROOM_WIDTH + 1, GRID_SPACING):
        for y in range(0, ROOM_HEIGHT + 1, GRID_SPACING):
            scanner_locations.append((x, y))

    # 3. Pre-calculate coverage for every possible placement
    print("Step 2: Pre-calculating coverage for every possible scanner placement...")
    coverage_map = {}
    point_to_idx = {point: i for i, point in enumerate(checkpoints)}

    for loc in scanner_locations:
        loc_x, loc_y = loc
        for name, specs in SCANNER_TYPES.items():
            covered_indices = set()
            if specs['shape'] == 'circle':
                radius_sq = specs['radius'] ** 2
                # Check points only within the circle's bounding box for efficiency
                min_x = max(0, math.floor(loc_x - specs['radius']))
                max_x = min(ROOM_WIDTH, math.ceil(loc_x + specs['radius']))
                min_y = max(0, math.floor(loc_y - specs['radius']))
                max_y = min(ROOM_HEIGHT, math.ceil(loc_y + specs['radius']))
                for px in range(min_x, max_x + 1):
                    for py in range(min_y, max_y + 1):
                        if (px - loc_x)**2 + (py - loc_y)**2 <= radius_sq:
                            covered_indices.add(point_to_idx.get((px, py)))
            elif specs['shape'] == 'square':
                half_side = specs['side'] / 2.0
                min_x = max(0, math.floor(loc_x - half_side))
                max_x = min(ROOM_WIDTH, math.ceil(loc_x + half_side))
                min_y = max(0, math.floor(loc_y - half_side))
                max_y = min(ROOM_HEIGHT, math.ceil(loc_y + half_side))
                for px in range(min_x, max_x + 1):
                     for py in range(min_y, max_y + 1):
                         covered_indices.add(point_to_idx.get((px, py)))
            covered_indices.discard(None) # Safely remove if a point was out of bounds
            coverage_map[(loc, name)] = covered_indices

    # 4. Greedy algorithm loop
    print("Step 3: Running greedy algorithm to find the most cost-effective scanner at each step...")
    placed_scanners = []
    total_cost = 0
    all_covered_points_indices = set()

    while len(all_covered_points_indices) < target_covered_points:
        best_cost_per_point = float('inf')
        best_placement = None
        best_newly_covered = set()

        for placement, covered_set in coverage_map.items():
            cost = SCANNER_TYPES[placement[1]]['cost']
            newly_covered_points = covered_set.difference(all_covered_points_indices)
            num_new = len(newly_covered_points)

            if num_new > 0:
                cost_per_point = cost / num_new
                if cost_per_point < best_cost_per_point:
                    best_cost_per_point = cost_per_point
                    best_placement = placement
                    best_newly_covered = newly_covered_points

        if best_placement is None:
            print("Warning: No more coverage can be added.")
            break

        # Place the best scanner found
        loc, name = best_placement
        cost = SCANNER_TYPES[name]['cost']

        placed_scanners.append({'type': name, 'location': loc, 'cost': cost})
        total_cost += cost
        all_covered_points_indices.update(best_newly_covered)

    # 5. Final Output
    print("\n--- Optimization Result ---")
    final_coverage = len(all_covered_points_indices) / total_points
    print(f"Final Coverage: {final_coverage:.4f} (Target was >= {COVERAGE_TARGET})")

    counts = {name: 0 for name in SCANNER_TYPES}
    cost_map = {name: specs['cost'] for name, specs in SCANNER_TYPES.items()}
    for scanner in placed_scanners:
        counts[scanner['type']] += 1

    print("\nFinal Scanner Configuration:")
    equation_parts_verbose = []
    equation_parts_numeric = []
    for name, count in sorted(counts.items()):
        if count > 0:
            cost = cost_map[name]
            sub_total = count * cost
            print(f"- Type {name}: {count} units")
            equation_parts_verbose.append(f"{count} * {cost}")
            equation_parts_numeric.append(str(sub_total))

    print("\nFinal Optimal Total Cost Calculation:")
    # Printing each number in the final equation
    print(f"Total Cost = {' + '.join(equation_parts_verbose)}")
    if len(equation_parts_numeric) > 1:
        print(f"           = {' + '.join(equation_parts_numeric)}")
    print(f"           = {total_cost}")

    print(f"\n<<<{total_cost}>>>")

solve_scanner_optimization()