import pulp
import math

def solve_museum_coverage():
    """
    Formulates and solves the museum scanner placement problem
    using Mixed-Integer Programming with the PuLP library.
    """
    # 1. Define constants
    ROOM_W = 140
    ROOM_H = 110
    COVERAGE_TARGET = 0.88
    GRID_STEP = 5

    scanners = {
        'C2': {'shape': 'circle', 'radius': 20, 'cost': 20000},
        'C1': {'shape': 'circle', 'radius': 5,  'cost': 1600}, # from 10m diameter
        'R1': {'shape': 'square', 'side': 10,   'cost': 2000},
    }
    scanner_types = list(scanners.keys())

    # 2. Set up the grid for placement and coverage evaluation
    grid_x = range(0, ROOM_W + 1, GRID_STEP)
    grid_y = range(0, ROOM_H + 1, GRID_STEP)
    grid_points = [(x, y) for x in grid_x for y in grid_y]

    total_points = len(grid_points)
    required_points_to_cover = math.ceil(total_points * COVERAGE_TARGET)

    print("--- Problem Setup ---")
    print(f"Room Dimensions: {ROOM_W}m x {ROOM_H}m")
    print(f"Grid Points (5mx5m cells): {total_points}")
    print(f"Coverage Target: {COVERAGE_TARGET*100}%, requiring at least {required_points_to_cover} points to be covered.\n")

    # 3. Pre-compute which placements cover which points
    point_covered_by = {p: [] for p in grid_points}
    for name, props in scanners.items():
        for sx, sy in grid_points:  # Scanner center location
            for px, py in grid_points:  # Target point location
                is_covered = False
                if props['shape'] == 'circle':
                    if (px - sx)**2 + (py - sy)**2 <= props['radius']**2:
                        is_covered = True
                elif props['shape'] == 'square':
                    half_side = props['side'] / 2
                    if abs(px - sx) <= half_side and abs(py - sy) <= half_side:
                        is_covered = True
                
                if is_covered:
                    point_covered_by[(px, py)].append((name, sx, sy))

    # 4. Formulate the MIP model using PuLP
    model = pulp.LpProblem("Museum_Scanner_Placement", pulp.LpMinimize)

    # Decision Variables
    scanner_vars = pulp.LpVariable.dicts(
        "Scanner",
        ((t, x, y) for t in scanner_types for (x, y) in grid_points),
        cat='Binary'
    )
    point_is_covered = pulp.LpVariable.dicts("IsCovered", grid_points, cat='Binary')

    # Objective Function: Minimize total cost
    cost_expr = pulp.lpSum(
        scanner_vars[t, loc[0], loc[1]] * scanners[t]['cost']
        for t in scanner_types for loc in grid_points
    )
    model += cost_expr, "Total_Cost"

    # Constraints
    # a. Total coverage constraint
    model += pulp.lpSum(point_is_covered.values()) >= required_points_to_cover, "TotalCoverage"

    # b. Link scanner placements to point coverage
    for p_loc, covering_scanners in point_covered_by.items():
        if covering_scanners:
             model += point_is_covered[p_loc] <= pulp.lpSum(
                [scanner_vars[s_type, sx, sy] for s_type, sx, sy in covering_scanners]
            ), f"Link_Coverage_at_{p_loc[0]}_{p_loc[1]}"
        else:
             model += point_is_covered[p_loc] == 0 # This point cannot be covered

    # c. At most one scanner per location
    for loc in grid_points:
        model += pulp.lpSum(scanner_vars[t, loc[0], loc[1]] for t in scanner_types) <= 1, f"OneScannerAt_{loc[0]}_{loc[1]}"

    # 5. Solve the problem
    print("--- Solving Optimization Problem ---")
    # PuLP will use its bundled CBC solver by default.
    model.solve()
    print("Solver finished.\n")

    # 6. Print the results
    print("--- Results ---")
    print(f"Optimization Status: {pulp.LpStatus[model.status]}")

    if pulp.LpStatus[model.status] == 'Optimal':
        total_cost = pulp.value(model.objective)
        print(f"\nOptimal Total Cost: {total_cost:.0f}")

        # Tally the results
        selected_scanners = {}
        for (t, x, y), var in scanner_vars.items():
            if var.varValue > 0.5:
                if t not in selected_scanners:
                    selected_scanners[t] = 0
                selected_scanners[t] += 1
        
        print("\nOptimal Scanner Configuration:")
        # Print in a consistent order
        for s_type in sorted(scanners.keys()):
            count = selected_scanners.get(s_type, 0)
            if count > 0:
                print(f"- Type {s_type}: {count} units")
        
        print("\nFinal Cost Equation:")
        cost_terms = []
        for s_type in sorted(scanners.keys()):
            count = selected_scanners.get(s_type, 0)
            if count > 0:
                cost_terms.append(f"{count} * {scanners[s_type]['cost']}")
        
        print(f"{' + '.join(cost_terms)} = {total_cost:.0f}")

        actual_covered_points = sum(p.varValue for p in point_is_covered.values())
        coverage_achieved = (actual_covered_points / total_points) * 100
        print(f"\nCoverage Achieved: {actual_covered_points:.0f}/{total_points} points ({coverage_achieved:.2f}%)")
    else:
        print("Could not find an optimal solution.")


if __name__ == '__main__':
    solve_museum_coverage()
