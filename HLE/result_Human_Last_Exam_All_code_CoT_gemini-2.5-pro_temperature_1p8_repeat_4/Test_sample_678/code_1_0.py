import pulp
import math
import sys

def solve_museum_coverage():
    """
    Formulates and solves the museum scanner coverage problem as a MILP.
    """
    # 1. Define constants and problem parameters
    ROOM_W = 140
    ROOM_H = 110
    GRID_STEP = 5

    TARGET_COVERAGE_RATIO = 0.88
    ROOM_AREA = ROOM_W * ROOM_H
    TARGET_AREA = ROOM_AREA * TARGET_COVERAGE_RATIO

    # Scanner types and their properties
    SCANNER_TYPES = {
        "C2": {"cost": 20000, "shape": "circle", "radius": 20, "side": None},
        "C1": {"cost": 1600,  "shape": "circle", "radius": 5, "side": None},
        "R1": {"cost": 2000,  "shape": "square", "side": 10, "radius": None}
    }
    
    # 2. Create a grid of points to model the area
    # Possible locations for scanners
    possible_locations = []
    for x in range(0, ROOM_W + 1, GRID_STEP):
        for y in range(0, ROOM_H + 1, GRID_STEP):
            possible_locations.append((x, y))

    # Grid points to be covered (same as locations in this model)
    points_to_cover = possible_locations

    print(f"Modeling the room with {len(points_to_cover)} discrete points...")

    # 3. Pre-calculate which scanners can cover which points
    # This avoids geometric calculations inside the optimization model
    point_coverage_map = {p: [] for p in points_to_cover}

    for l_x, l_y in possible_locations:
        loc = (l_x, l_y)
        for name, props in SCANNER_TYPES.items():
            if props["shape"] == "circle":
                radius_sq = props["radius"]**2
                for p_x, p_y in points_to_cover:
                    if (p_x - l_x)**2 + (p_y - l_y)**2 <= radius_sq:
                        point_coverage_map[(p_x, p_y)].append((name, loc))
            elif props["shape"] == "square":
                half_side = props["side"] / 2
                for p_x, p_y in points_to_cover:
                    if abs(p_x - l_x) <= half_side and abs(p_y - l_y) <= half_side:
                        point_coverage_map[(p_x, p_y)].append((name, loc))

    # 4. Set up the MILP problem
    prob = pulp.LpProblem("Museum_Coverage_Optimization", pulp.LpMinimize)

    # Decision Variables
    # x_loc_type = 1 if scanner of 'type' is placed at 'loc'
    scanner_vars = pulp.LpVariable.dicts(
        "Scanner",
        ((loc, type_name) for loc in possible_locations for type_name in SCANNER_TYPES),
        cat='Binary'
    )
    # y_point = 1 if 'point' is covered
    covered_point_vars = pulp.LpVariable.dicts(
        "PointCovered",
        points_to_cover,
        cat='Binary'
    )

    # 5. Define the Objective Function (minimize total cost)
    total_cost = pulp.lpSum(
        scanner_vars[loc, type_name] * SCANNER_TYPES[type_name]["cost"]
        for loc in possible_locations for type_name in SCANNER_TYPES
    )
    prob += total_cost

    # 6. Define Constraints
    # a) Total Coverage Constraint: at least 88% area coverage
    area_per_point = GRID_STEP * GRID_STEP
    num_points_needed = math.ceil(TARGET_AREA / area_per_point)
    prob += pulp.lpSum(covered_point_vars[p] for p in points_to_cover) >= num_points_needed
    
    # b) Link scanner placement to point coverage
    for point, covering_scanners in point_coverage_map.items():
        prob += pulp.lpSum(scanner_vars[loc, type_name] for type_name, loc in covering_scanners) >= covered_point_vars[point]

    # c) At most one scanner per location
    for loc in possible_locations:
        prob += pulp.lpSum(scanner_vars[loc, type_name] for type_name in SCANNER_TYPES) <= 1

    # 7. Solve the problem
    print("Solving the optimization problem. This may take a few minutes...")
    # Using the default CBC solver that comes with pulp
    prob.solve()

    # 8. Output the results
    print("\n--- Optimization Results ---")
    status = pulp.LpStatus[prob.status]
    print(f"Status: {status}")

    if status == 'Optimal':
        final_cost = pulp.value(prob.objective)
        
        counts = {name: 0 for name in SCANNER_TYPES}
        for var in scanner_vars.values():
            if var.varValue > 0.5:
                loc, type_name = var.name.split('_(')[1].strip(')').replace("'", "").replace(" ", "").split(',')
                counts[loc] += 1
        
        # Build and print the final cost equation
        cost_breakdown = []
        for name, num in counts.items():
            if num > 0:
                cost = SCANNER_TYPES[name]["cost"]
                cost_breakdown.append(f"{num} * {cost}")
        
        print("\nOptimal Scanner Configuration:")
        equation_str = " + ".join(cost_breakdown)
        print(f"Cost Equation: {equation_str} = {int(final_cost)}")

        # Print details
        total_points_covered = sum(p.varValue for p in covered_point_vars.values())
        actual_coverage_area = total_points_covered * area_per_point
        actual_coverage_ratio = actual_coverage_area / ROOM_AREA
        
        print(f"Optimal Total Cost: {int(final_cost)}")
        print(f"Coverage Achieved: {actual_coverage_area:.2f} / {ROOM_AREA} sq m ({actual_coverage_ratio:.2%})")
        print("Number of each scanner type used:")
        for name, num in counts.items():
            print(f"- {name}: {num}")

        # Final answer format
        print(f"\n<<<${int(final_cost)}>>>")
        
    else:
        print("Optimal solution not found.")


if __name__ == '__main__':
    try:
        import pulp
    except ImportError:
        print("Error: The 'pulp' library is required. Please install it using 'pip install pulp'")
        sys.exit(1)
    solve_museum_coverage()