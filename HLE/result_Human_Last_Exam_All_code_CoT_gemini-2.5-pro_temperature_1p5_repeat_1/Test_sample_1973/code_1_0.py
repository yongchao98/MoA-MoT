import pulp

def solve_unicorn_problem():
    """
    Finds the minimum number of unicorns to attack all black squares on an 8x8x8 board
    using Integer Linear Programming.
    """
    N = 8

    # 1. Create the ILP problem
    model = pulp.LpProblem("Unicorn_Cover_Problem", pulp.LpMinimize)

    # 2. Define Decision Variables
    # x_i_j_k = 1 if a unicorn is placed on square (i,j,k), 0 otherwise
    square_indices = [(i, j, k) for i in range(N) for j in range(N) for k in range(N)]
    unicorn_vars = pulp.LpVariable.dicts("Unicorn", square_indices, cat='Binary')

    # 3. Define the Objective Function
    # Minimize the total number of unicorns placed
    model += pulp.lpSum(unicorn_vars[i] for i in square_indices)

    # 4. Define the Constraints
    # Every black square must be attacked by at least one unicorn.
    black_squares = []
    for r in range(N):
        for c in range(N):
            for l in range(N):
                # A square is black if the sum of its coordinates is odd.
                if (r + c + l) % 2 != 0:
                    black_squares.append((r, c, l))

    # For each black square, add a constraint
    for r_b, c_b, l_b in black_squares:
        attacking_unicorns = []
        # Find all squares (i,j,k) from which a unicorn can attack the black square
        for r_u in range(N):
            for c_u in range(N):
                for l_u in range(N):
                    # A unicorn at (r_u, c_u, l_u) attacks (r_b, c_b, l_b) if they
                    # are on a 3D diagonal.
                    if abs(r_b - r_u) == abs(c_b - c_u) == abs(l_b - l_u):
                        attacking_unicorns.append(unicorn_vars[(r_u, c_u, l_u)])
        
        # The sum of unicorns that can attack this black square must be >= 1
        model += pulp.lpSum(attacking_unicorns) >= 1, f"Attack_Constraint_({r_b},{c_b},{l_b})"

    # 5. Solve the problem
    # The default CBC solver that comes with pulp is sufficient for this.
    model.solve()

    # 6. Print the results
    print(f"Solver Status: {pulp.LpStatus[model.status]}")
    
    if model.status == pulp.LpStatusOptimal:
        min_unicorns = int(pulp.value(model.objective))
        print(f"\nMinimum number of unicorns needed: {min_unicorns}")
        
        print("\nThe final equation is constructed by placing one unicorn at each of the following locations:")
        
        unicorn_placements = []
        for i in square_indices:
            if unicorn_vars[i].varValue == 1:
                unicorn_placements.append(str(i))
                print(f"  {i}")
        
        equation_str = " + ".join(["1"] * len(unicorn_placements))
        print(f"\nEquation: {equation_str} = {min_unicorns}")
        
        # This is the final numerical answer as requested by the user prompt format.
        print(f"\n<<<{min_unicorns}>>>")

if __name__ == "__main__":
    solve_unicorn_problem()