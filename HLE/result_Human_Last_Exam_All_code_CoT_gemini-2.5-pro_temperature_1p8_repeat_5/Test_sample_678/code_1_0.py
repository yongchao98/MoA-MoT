import numpy as np
from scipy.optimize import milp
import time

def solve_scanner_placement():
    """
    Formulates and solves the museum scanner placement problem as a MILP.
    """
    # 1. Constants and Setup
    ROOM_W, ROOM_H = 140, 110
    TARGET_COVERAGE_RATIO = 0.88
    GRID_STEP = 5

    SCANNER_SPECS = {
        'C2': {'type': 'circle', 'param': 20, 'cost': 20000}, # radius 20m
        'C1': {'type': 'circle', 'param': 5, 'cost': 1600},   # diameter 10m -> radius 5m
        'R1': {'type': 'square', 'param': 5, 'cost': 2000},   # side 10m -> half-side 5m
    }
    SCANNER_TYPES = list(SCANNER_SPECS.keys())
    N_TYPES = len(SCANNER_TYPES)

    # 2. Grid Setup
    x_coords = np.arange(0, ROOM_W + 1, GRID_STEP)
    y_coords = np.arange(0, ROOM_H + 1, GRID_STEP)
    
    locations = []
    for x in x_coords:
        for y in y_coords:
            locations.append((x, y))
            
    N_LOCATIONS = len(locations)
    # Use the same grid for placement locations and coverage test points
    test_points = locations
    N_POINTS = len(test_points)

    TARGET_COVERED_POINTS = int(np.ceil(TARGET_COVERAGE_RATIO * N_POINTS))
    
    print(f"Room size: {ROOM_W}x{ROOM_H}m. Placement grid: {len(x_coords)}x{len(y_coords)} ({N_LOCATIONS} locations).")
    print(f"Coverage requirement: {TARGET_COVERAGE_RATIO:.0%} of {N_POINTS} test points = {TARGET_COVERED_POINTS} points.\n")


    # 3. Pre-computation of Coverage Matrix
    # covers_matrix[type_idx, loc_idx, point_idx] is True if scanner of type at loc covers point
    print("Pre-calculating coverage matrix for all scanner types and locations...")
    covers_matrix = np.zeros((N_TYPES, N_LOCATIONS, N_POINTS), dtype=bool)

    for type_idx, type_name in enumerate(SCANNER_TYPES):
        scanner = SCANNER_SPECS[type_name]
        for loc_idx, (lx, ly) in enumerate(locations):
            for pt_idx, (px, py) in enumerate(test_points):
                if scanner['type'] == 'circle':
                    radius_sq = scanner['param']**2
                    if (lx - px)**2 + (ly - py)**2 <= radius_sq:
                        covers_matrix[type_idx, loc_idx, pt_idx] = True
                elif scanner['type'] == 'square':
                    half_side = scanner['param']
                    if abs(lx - px) <= half_side and abs(ly - py) <= half_side:
                        covers_matrix[type_idx, loc_idx, pt_idx] = True
    print("Coverage matrix calculation complete.\n")


    # 4. MILP Formulation
    # Variable vector 'x' structure:
    # [s_C2_0, ..., s_C2_N-1,  (N_LOCATIONS vars for C2 placement)
    #  s_C1_0, ..., s_C1_N-1,  (N_LOCATIONS vars for C1 placement)
    #  s_R1_0, ..., s_R1_N-1,  (N_LOCATIONS vars for R1 placement)
    #  p_0,    ..., p_M-1]     (N_POINTS vars for point coverage status)
    
    N_SCANNER_VARS = N_TYPES * N_LOCATIONS
    N_TOTAL_VARS = N_SCANNER_VARS + N_POINTS

    # Objective function (c): Minimize total cost
    c = np.zeros(N_TOTAL_VARS)
    for type_idx, type_name in enumerate(SCANNER_TYPES):
        cost = SCANNER_SPECS[type_name]['cost']
        start_idx = type_idx * N_LOCATIONS
        end_idx = start_idx + N_LOCATIONS
        c[start_idx:end_idx] = cost
    # The coverage variables p_j have zero cost in the objective function

    # Constraints (in the form A_ub @ x <= b_ub)
    A_ub_rows = []
    b_ub_vals = []

    # Constraint 1: At most one scanner per location i
    # s_C2_i + s_C1_i + s_R1_i <= 1
    for i in range(N_LOCATIONS):
        row = np.zeros(N_TOTAL_VARS)
        row[i] = 1                              # C2 at loc i
        row[i + N_LOCATIONS] = 1                # C1 at loc i
        row[i + 2 * N_LOCATIONS] = 1            # R1 at loc i
        A_ub_rows.append(row)
        b_ub_vals.append(1)

    # Constraint 2: Link scanner placement to point coverage
    # p_j <= sum over all scanners k,i (s_k_i * covers_matrix[k, i, j])
    # Rewritten as: p_j - sum(s_k_i * covers_matrix[k,i,j]) <= 0
    for j in range(N_POINTS): # for each test point j
        row = np.zeros(N_TOTAL_VARS)
        row[N_SCANNER_VARS + j] = 1 # Coeff for p_j is 1
        for k in range(N_TYPES): # for each scanner type
            for i in range(N_LOCATIONS): # for each location
                if covers_matrix[k, i, j]:
                    scanner_var_idx = k * N_LOCATIONS + i
                    row[scanner_var_idx] = -1 # Coeff for s_k_i is -1
        A_ub_rows.append(row)
        b_ub_vals.append(0)

    # Constraint 3: Total coverage requirement
    # sum(p_j for all j) >= TARGET_COVERED_POINTS
    # Rewritten as: -sum(p_j) <= -TARGET_COVERED_POINTS
    row = np.zeros(N_TOTAL_VARS)
    row[N_SCANNER_VARS:] = -1 # -1 for all p_j variables
    A_ub_rows.append(row)
    b_ub_vals.append(-TARGET_COVERED_POINTS)

    A_ub = np.array(A_ub_rows)
    b_ub = np.array(b_ub_vals)

    # All variables are binary (integer and between 0 and 1)
    integrality = np.ones(N_TOTAL_VARS)
    bounds = (0, 1)

    # 5. Solve the MILP
    print("Solving the optimization problem. This may take a few minutes...")
    start_time = time.time()
    res = milp(c=c, integrality=integrality, bounds=bounds, constraints={'A_ub': A_ub, 'b_ub': b_ub})
    end_time = time.time()
    print(f"Solver finished in {end_time - start_time:.2f} seconds.\n")

    # 6. Output the Results
    if res.success:
        total_cost = res.fun
        solution_vars = np.round(res.x).astype(int)
        
        counts = {}
        for type_idx, type_name in enumerate(SCANNER_TYPES):
            start_idx = type_idx * N_LOCATIONS
            end_idx = start_idx + N_LOCATIONS
            counts[type_name] = np.sum(solution_vars[start_idx:end_idx])

        covered_points = np.sum(solution_vars[N_SCANNER_VARS:])
        coverage_ratio = covered_points / N_POINTS

        print("--- Optimal Solution Found ---")
        print(f"Minimum Total Cost: {total_cost:.0f}")
        print(f"Resulting Coverage: {covered_points}/{N_POINTS} points ({coverage_ratio:.2%}), meets >= {TARGET_COVERAGE_RATIO:.0%} target.")
        print("\nScanner quantities and costs:")
        
        cost_equation_parts = []
        for type_name, count in counts.items():
            if count > 0:
                cost = SCANNER_SPECS[type_name]['cost']
                subtotal = count * cost
                print(f"  - Type {type_name}: {count} scanners x {cost} = {subtotal}")
                cost_equation_parts.append(f"{count} * {cost}")

        equation_str = " + ".join(cost_equation_parts)
        print(f"\nFinal Cost Equation:")
        print(f"{equation_str} = {int(total_cost)}")

    else:
        print("Solver did not find an optimal solution.")
        print(f"Solver status: {res.message}")
        total_cost = -1 # Indicate failure

    return int(total_cost)

if __name__ == '__main__':
    final_cost = solve_scanner_placement()
    if final_cost != -1:
        print(f"\n<<<${final_cost}>>>")
