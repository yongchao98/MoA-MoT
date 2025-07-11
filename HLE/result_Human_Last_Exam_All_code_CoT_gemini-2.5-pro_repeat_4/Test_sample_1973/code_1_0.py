import pulp

def solve_unicorn_problem():
    """
    Calculates the minimum number of unicorns to attack all black squares on an 8x8x8 board
    using Integer Linear Programming.
    """
    N = 8

    # The problem decomposes into 4 identical subproblems. We solve one:
    # covering OOO (Odd,Odd,Odd) squares with unicorns on EEE (Even,Even,Even) squares.

    # 1. Define the sets of squares for the subproblem
    unicorn_positions = []
    for x in range(N):
        for y in range(N):
            for z in range(N):
                if x % 2 == 0 and y % 2 == 0 and z % 2 == 0:
                    unicorn_positions.append((x, y, z))

    target_squares = []
    for x in range(N):
        for y in range(N):
            for z in range(N):
                if x % 2 != 0 and y % 2 != 0 and z % 2 != 0:
                    target_squares.append((x, y, z))

    # 2. Set up the ILP problem
    prob = pulp.LpProblem("Unicorn_Cover_Subproblem", pulp.LpMinimize)

    # 3. Define decision variables
    # A dictionary of binary variables, one for each possible unicorn position
    unicorn_vars = pulp.LpVariable.dicts("Unicorn", unicorn_positions, cat='Binary')

    # 4. Define the objective function
    # Minimize the total number of unicorns placed
    prob += pulp.lpSum(unicorn_vars[p] for p in unicorn_positions)

    # 5. Define the constraints
    # For each target square, at least one unicorn must be able to attack it.
    for t in target_squares:
        attacking_unicorns = []
        for p in unicorn_positions:
            dx = abs(p[0] - t[0])
            dy = abs(p[1] - t[1])
            dz = abs(p[2] - t[2])
            # Check for unicorn attack condition
            if dx > 0 and dx == dy and dx == dz:
                attacking_unicorns.append(unicorn_vars[p])
        
        # Add the constraint for this target square
        prob += pulp.lpSum(attacking_unicorns) >= 1, f"Cover_{t}"

    # 6. Solve the ILP problem
    # The solver will run quietly (msg=0)
    solver = pulp.PULP_CBC_CMD(msg=0)
    prob.solve(solver)

    # 7. Process and print the results
    print(f"Solving one of the four symmetric subproblems...")
    print("-" * 50)
    
    if pulp.LpStatus[prob.status] == 'Optimal':
        subproblem_solution = int(pulp.value(prob.objective))
        
        print(f"Minimum unicorns for one subproblem: {subproblem_solution}")
        print("Unicorns are placed at the following EEE coordinates:")
        
        chosen_unicorns = []
        for p in unicorn_positions:
            if unicorn_vars[p].varValue > 0.5:
                chosen_unicorns.append(p)
                print(p)
        
        # Print the equation for the subproblem
        ones = ["1"] * len(chosen_unicorns)
        equation = " + ".join(ones)
        print(f"\nThe equation for this subproblem is: {equation} = {subproblem_solution}")
        
        print("-" * 50)
        # Calculate and print the final result for the full 8x8x8 board
        total_unicorns = 4 * subproblem_solution
        print("The full problem consists of 4 such independent subproblems.")
        final_equation = f"{subproblem_solution} + {subproblem_solution} + {subproblem_solution} + {subproblem_solution} = {total_unicorns}"
        print(f"The final equation is: {final_equation}")
        print(f"\nThe minimum number of unicorns to attack all black squares is: {total_unicorns}")

    else:
        print(f"Could not find an optimal solution. Status: {pulp.LpStatus[prob.status]}")

if __name__ == '__main__':
    solve_unicorn_problem()
