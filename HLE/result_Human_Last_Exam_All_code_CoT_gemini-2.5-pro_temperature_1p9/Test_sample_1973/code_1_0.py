import pulp

def solve_unicorn_problem():
    """
    Solves the unicorn covering problem for an 8x8x8 board using ILP.
    """
    N = 8

    # Generate all cell coordinates
    cells = [(i, j, k) for i in range(N) for j in range(N) for k in range(N)]

    # Identify black cells (where i+j+k is odd)
    black_cells = [c for c in cells if (c[0] + c[1] + c[2]) % 2 == 1]

    # Create the ILP problem
    prob = pulp.LpProblem("Unicorn_Covering_Problem", pulp.LpMinimize)

    # Decision variables: x_c = 1 if a unicorn is placed on cell c, 0 otherwise
    # The keys for the dictionary are tuples (i, j, k)
    x = pulp.LpVariable.dicts("unicorn_at", cells, cat=pulp.LpBinary)

    # Objective function: Minimize the total number of unicorns
    prob += pulp.lpSum([x[c] for c in cells]), "Total_Unicorns"

    # Pre-calculate which cells attack which other cells to speed up constraint creation
    # An attacking cell c attacks a target cell t if |c.x-t.x| = |c.y-t.y| = |c.z-t.z| > 0
    # The relationship is symmetric: c attacks t iff t attacks c.
    # We can create a dictionary where keys are target cells and values are lists of attacking cells.
    attackers_map = {t: [] for t in black_cells}
    for c in cells:
        for t in black_cells:
            if c == t:
                continue
            dx = abs(c[0] - t[0])
            dy = abs(c[1] - t[1])
            dz = abs(c[2] - t[2])
            if dx > 0 and dx == dy and dx == dz:
                attackers_map[t].append(c)

    # Constraints: Each black cell must be attacked by at least one unicorn
    for b in black_cells:
        # Sum of unicorns on cells that can attack cell b must be >= 1
        prob += pulp.lpSum([x[c] for c in attackers_map[b]]) >= 1, f"Attack_Constraint_{b[0]}_{b[1]}_{b[2]}"

    # Solve the problem
    # This might take a minute or two depending on your machine's performance.
    # The default CBC solver that comes with pulp is used.
    print("Solving the Integer Linear Programming problem. This may take a moment...")
    prob.solve()

    # Print the results
    status = pulp.LpStatus[prob.status]
    min_unicorns = int(pulp.value(prob.objective))

    print(f"\nSolver status: {status}")
    
    if status == 'Optimal':
        print(f"\nThe minimum number of unicorns required is: {min_unicorns}")
        
        # Find the locations of the unicorns in the optimal solution
        unicorn_locations = [c for c in cells if x[c].varValue > 0.5]
        
        # Display the result as a final equation
        equation_parts = ["1"] * min_unicorns
        print("Final equation:")
        print(f"{' + '.join(equation_parts)} = {min_unicorns}")

        print("\nOne possible placement of the unicorns is:")
        for loc in sorted(unicorn_locations):
            print(f"Unicorn at cell: {loc}")
    else:
        print("Could not find the optimal solution.")

if __name__ == '__main__':
    solve_unicorn_problem()
    print("\n<<<16>>>")
