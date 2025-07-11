import pulp
import math

def solve_museum_coverage():
    """
    Formulates and solves the museum scanner placement problem as a Mixed-Integer Program.
    """
    # 1. Define the Problem Space
    ROOM_WIDTH = 140
    ROOM_HEIGHT = 110
    GRID_STEP = 5
    TARGET_COVERAGE = 0.88

    # Scanner definitions: type -> {property: value}
    SCANNERS = {
        'C2': {'shape': 'circle', 'range': 20, 'cost': 20000}, # range is radius
        'C1': {'shape': 'circle', 'range': 5,  'cost': 1600},  # d=10m -> r=5m
        'R1': {'shape': 'square', 'range': 5,  'cost': 2000},  # side=10m -> range from center is 5m
    }

    # Generate the grid of possible locations (and points to be covered)
    x_coords = range(0, ROOM_WIDTH + 1, GRID_STEP)
    y_coords = range(0, ROOM_HEIGHT + 1, GRID_STEP)
    locations = [(x, y) for x in x_coords for y in y_coords]
    num_locations = len(locations)
    target_covered_points = math.floor(num_locations * TARGET_COVERAGE)

    # 2. Pre-computation: Map points to the placements that can cover them
    print("Pre-computing coverage map... (this may take a moment)")
    point_to_placements = {p: [] for p in locations}
    for s_type, props in SCANNERS.items():
        for l in locations:
            lx, ly = l
            for p in locations:
                px, py = p
                is_covered = False
                if props['shape'] == 'circle':
                    # Use squared distance to avoid sqrt
                    if (px - lx)**2 + (py - ly)**2 <= props['range']**2:
                        is_covered = True
                elif props['shape'] == 'square':
                    if abs(px - lx) <= props['range'] and abs(py - ly) <= props['range']:
                        is_covered = True
                
                if is_covered:
                    point_to_placements[p].append((s_type, l))
    print("Pre-computation finished.")

    # 3. Formulate the MIP problem
    prob = pulp.LpProblem("Museum_Scanner_Placement", pulp.LpMinimize)

    # Decision Variables
    # placement_vars[s_type, location] = 1 if scanner s_type is at location, 0 otherwise
    placement_vars = pulp.LpVariable.dicts("placement",
                                           ((s_type, loc) for s_type in SCANNERS for loc in locations),
                                           cat='Binary')
    
    # covered_vars[location] = 1 if the point at location is covered, 0 otherwise
    covered_vars = pulp.LpVariable.dicts("covered", locations, cat='Binary')

    # 4. Define the Objective Function (minimize total cost)
    total_cost = pulp.lpSum(SCANNERS[s_type]['cost'] * placement_vars[(s_type, loc)]
                             for s_type in SCANNERS for loc in locations)
    prob += total_cost

    # 5. Define the Constraints
    # Constraint 1: Achieve target coverage
    prob += pulp.lpSum(covered_vars[loc] for loc in locations) >= target_covered_points, "TotalCoverage"

    # Constraint 2: Link coverage to placements
    # A point is covered only if at least one active scanner covers it.
    for p in locations:
        prob += pulp.lpSum(placement_vars[(s_type, l)] for (s_type, l) in point_to_placements[p]) >= covered_vars[p], f"LinkCoverage_{p}"

    # Constraint 3: At most one scanner per location
    for loc in locations:
        prob += pulp.lpSum(placement_vars[(s_type, loc)] for s_type in SCANNERS) <= 1, f"OneScannerPerLocation_{loc}"
    
    # 6. Solve the problem
    print("Solving the optimization problem...")
    # Using the default CBC solver. msg=1 shows solver output.
    solver = pulp.PULP_CBC_CMD(msg=0)
    prob.solve(solver)
    print("Solver finished.")

    # 7. Process and print the results
    if pulp.LpStatus[prob.status] == 'Optimal':
        final_cost = pulp.value(prob.objective)
        scanner_counts = {s_type: 0 for s_type in SCANNERS}
        for s_type in SCANNERS:
            for loc in locations:
                if placement_vars[(s_type, loc)].varValue > 0.5:
                    scanner_counts[s_type] += 1
        
        c2_count = scanner_counts.get('C2', 0)
        c1_count = scanner_counts.get('C1', 0)
        r1_count = scanner_counts.get('R1', 0)

        c2_cost = SCANNERS['C2']['cost']
        c1_cost = SCANNERS['C1']['cost']
        r1_cost = SCANNERS['R1']['cost']
        
        print("\nOptimal Solution Found:")
        print(f"The minimum cost to achieve at least {TARGET_COVERAGE*100}% coverage is: {int(final_cost)}")
        print("\nThis cost is achieved using:")
        print(f" - {c2_count} of C2 Scanners")
        print(f" - {c1_count} of C1 Scanners")
        print(f" - {r1_count} of R1 Scanners")

        print("\nFinal Cost Calculation:")
        print(f"{c2_count} * {c2_cost} + {c1_count} * {c1_cost} + {r1_count} * {r1_cost} = {int(final_cost)}")
        
        print(f"\n<<<{int(final_cost)}>>>")

    else:
        print(f"Could not find an optimal solution. Status: {pulp.LpStatus[prob.status]}")

if __name__ == '__main__':
    solve_museum_coverage()