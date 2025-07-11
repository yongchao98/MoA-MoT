import pulp

def solve_unicorn_problem():
    """
    Calculates the minimum number of unicorns to attack all black squares on an 8x8x8 board
    by solving a decomposed subproblem using Integer Linear Programming.
    """
    N = 8

    # --- Step 1: Define the subproblem ---
    # The problem decomposes into 4 independent sets based on parity invariants.
    # We solve for one set: (x+y)%2 == 0 and (y+z)%2 == 0.
    # This implies x, y, and z all have the same parity.
    
    # Potential unicorn positions are all cells in this set.
    subproblem_cells = []
    for x in range(N):
        for y in range(N):
            for z in range(N):
                if (x % 2 == y % 2) and (y % 2 == z % 2):
                    subproblem_cells.append((x, y, z))

    # The cells to be covered are the black cells within this set.
    # A cell (x,y,z) is black if (x+y+z) is odd.
    # If x,y,z have the same parity, their sum is odd only if all three are odd.
    black_cells_to_cover = [c for c in subproblem_cells if (c[0] + c[1] + c[2]) % 2 == 1]

    # --- Step 2: Pre-calculate attack patterns ---
    # For each black cell q, find all cells p in our set that can attack it.
    attacks_map = {}
    for q_cell in black_cells_to_cover:
        qx, qy, qz = q_cell
        attackers = []
        for p_cell in subproblem_cells:
            px, py, pz = p_cell
            
            # A unicorn at a cell attacks that cell itself.
            if p_cell == q_cell:
                attackers.append(p_cell)
                continue
            
            # A unicorn attacks along diagonals.
            dx = abs(qx - px)
            dy = abs(qy - py)
            dz = abs(qz - pz)

            if dx > 0 and dx == dy and dx == dz:
                attackers.append(p_cell)
        
        attacks_map[q_cell] = attackers

    # --- Step 3: Formulate and solve the ILP problem ---
    # Create the ILP model
    prob = pulp.LpProblem("Unicorn_Subproblem", pulp.LpMinimize)

    # Create binary decision variables for each potential unicorn position
    unicorn_vars = pulp.LpVariable.dicts("Unicorn", subproblem_cells, cat='Binary')

    # Set the objective function: minimize the number of unicorns
    prob += pulp.lpSum(unicorn_vars[c] for c in subproblem_cells)

    # Add constraints: each black cell must be attacked by at least one unicorn
    for q_cell in black_cells_to_cover:
        prob += pulp.lpSum(unicorn_vars[p] for p in attacks_map[q_cell]) >= 1, f"Cover_{q_cell}"

    # Solve the problem (suppress solver output)
    prob.solve(pulp.PULP_CBC_CMD(msg=False))

    # --- Step 4: Output the final result ---
    subproblem_solution = int(pulp.value(prob.objective))
    total_unicorns = 4 * subproblem_solution

    print(f"The 8x8x8 board decomposes into 4 independent subproblems.")
    print(f"The minimum number of unicorns required to solve one subproblem is: {subproblem_solution}")
    print(f"The final equation is: 4 * {subproblem_solution} = {total_unicorns}")
    print(f"Therefore, the minimum number of unicorns to attack all black squares is {total_unicorns}.")
    
    return total_unicorns

if __name__ == '__main__':
    result = solve_unicorn_problem()
    print(f"\n<<<solution>>>\n{result}")
