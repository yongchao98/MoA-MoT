import pulp

def solve_unicorn_problem():
    """
    Calculates the minimum number of unicorns to attack all black squares
    on an 8x8x8 chessboard using Integer Linear Programming.
    """
    N = 8
    problem_name = "Unicorn_Dominating_Set"

    # Step 1: Identify all black squares on the 8x8x8 board.
    # A square (x, y, z) is black if (x + y + z) is odd.
    black_squares = []
    for x in range(N):
        for y in range(N):
            for z in range(N):
                if (x + y + z) % 2 == 1:
                    black_squares.append((x, y, z))

    # Step 2: Create the ILP problem.
    # We want to minimize the number of unicorns.
    prob = pulp.LpProblem(problem_name, pulp.LpMinimize)

    # Step 3: Create the decision variables.
    # We only need to consider placing unicorns on black squares, as a unicorn
    # on a white square cannot attack a black square.
    # unicorn_vars[(x, y, z)] = 1 if a unicorn is on (x,y,z), 0 otherwise.
    unicorn_vars = pulp.LpVariable.dicts(
        "Unicorn", black_squares, cat=pulp.LpBinary
    )

    # Step 4: Set the objective function.
    # Minimize the sum of all unicorn variables.
    prob += pulp.lpSum(unicorn_vars[s] for s in black_squares), "Total_Unicorns"

    # Helper function to determine if a unicorn at s1 can attack s2.
    def can_attack(s1, s2):
        """A unicorn moves along 3D diagonals."""
        return abs(s1[0] - s2[0]) == abs(s1[1] - s2[1]) == abs(s1[2] - s2[2])

    # Step 5: Add the constraints.
    # For each black square, it must be attacked by at least one unicorn.
    for target_square in black_squares:
        # The set of attackers includes any unicorn that can move to the target square.
        # This also includes a unicorn placed on the target square itself.
        attackers = [
            unicorn_vars[source_square]
            for source_square in black_squares
            if can_attack(source_square, target_square)
        ]
        # The sum of attacking unicorns must be at least 1.
        prob += pulp.lpSum(attackers) >= 1, f"Attack_Constraint_{target_square}"

    # Step 6: Solve the ILP problem.
    # PuLP uses the CBC solver by default.
    prob.solve()

    # Step 7: Print the results.
    print(f"Solver status: {pulp.LpStatus[prob.status]}")
    if prob.status == pulp.LpStatusOptimal:
        min_unicorns = int(pulp.value(prob.objective))
        print(f"The minimum number of unicorns needed is: {min_unicorns}")
        
        # As requested, create and print the final equation.
        equation_parts = ["1"] * min_unicorns
        equation_str = " + ".join(equation_parts)
        print(f"The final equation is: {equation_str} = {min_unicorns}")

        # The following lines can be uncommented to see the exact positions
        # of the unicorns in one possible optimal solution.
        # print("\nOne possible placement for the unicorns:")
        # for s in black_squares:
        #     if unicorn_vars[s].varValue == 1:
        #         print(f"  - {s}")
    else:
        print("Could not find the optimal solution.")


if __name__ == "__main__":
    solve_unicorn_problem()