import pulp

def solve_unicorn_domination():
    """
    Calculates the minimum number of unicorns to attack all black squares
    on an 8x8x8 3D chessboard using Integer Linear Programming.
    """
    N = 8
    problem_name = "Unicorn_Dominating_Set"

    # 1. Initialize the ILP model
    # We want to minimize the number of unicorns.
    model = pulp.LpProblem(problem_name, pulp.LpMinimize)

    # 2. Define the board and identify black squares
    all_cells = [(i, j, k) for i in range(N) for j in range(N) for k in range(N)]
    
    # A cell (i,j,k) is black if the sum of its coordinates is odd.
    black_cells = [cell for cell in all_cells if (cell[0] + cell[1] + cell[2]) % 2 != 0]

    # 3. Define the decision variables
    # We create a binary variable for each cell on the board.
    # x_cell = 1 if a unicorn is placed at 'cell', 0 otherwise.
    x = pulp.LpVariable.dicts("unicorn_at", all_cells, cat=pulp.LpBinary)

    # 4. Set the objective function
    # The objective is to minimize the total number of unicorns placed.
    model += pulp.lpSum(x[cell] for cell in all_cells), "Total_Unicorns"

    # 5. Add the constraints
    # For each black cell, we must ensure it is attacked by at least one unicorn.
    for b_cell in black_cells:
        r, c, l = b_cell
        
        # An attacking unicorn can be on any cell 'p_cell' (including b_cell itself)
        # that lies on a space diagonal with 'b_cell'.
        # This condition is |r-i| == |c-j| == |l-k|.
        
        attacking_cells = []
        for p_cell in all_cells:
            i, j, k = p_cell
            if abs(r - i) == abs(c - j) == abs(l - k):
                attacking_cells.append(p_cell)
        
        # The constraint: The sum of unicorns on attacking cells must be >= 1.
        model += pulp.lpSum(x[p] for p in attacking_cells) >= 1, f"attack_coverage_for_{b_cell}"

    # 6. Solve the model
    # The solver will find the minimum value for the objective function.
    # PuLP uses the CBC solver by default.
    model.solve()

    # 7. Print the final result
    # The final equation would be the sum of the unicorn variables that equal 1.
    # For example, if the answer is 16, the equation is 1 + 1 + ... (16 times) = 16.
    # We will print the final numerical result.
    min_unicorns = int(pulp.value(model.objective))
    
    print("Raumschach (8x8x8) Unicorn Domination Problem")
    print(f"Status: {pulp.LpStatus[model.status]}")
    print("\nTo attack all 256 black squares, the final equation for the minimum number of unicorns is:")
    
    # Create the string for the final "equation"
    equation_str = " + ".join(["1"] * min_unicorns)
    print(f"{equation_str} = {min_unicorns}")
    
    print(f"\nThe minimum number of unicorns required is: {min_unicorns}")

if __name__ == '__main__':
    solve_unicorn_domination()
