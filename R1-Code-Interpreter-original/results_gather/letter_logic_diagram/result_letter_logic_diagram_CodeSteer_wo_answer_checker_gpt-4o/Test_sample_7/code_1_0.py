from constraint import Problem, AllDifferentConstraint

def solve_puzzle():
    problem = Problem()

    # Define variables for each cell in the grid
    variables = [(r, c) for r in range(7) for c in range(7)]
    for var in variables:
        problem.addVariable(var, 'abcdefg')

    # Add constraints for each row
    for r in range(7):
        problem.addConstraint(AllDifferentConstraint(), [(r, c) for c in range(7)])

    # Add constraints for each column
    for c in range(7):
        problem.addConstraint(AllDifferentConstraint(), [(r, c) for r in range(7)])

    # Add constraint for the minor diagonal
    minor_diagonal = [(i, 6-i) for i in range(7)]
    problem.addConstraint(lambda *args: len(set(args)) == 1, minor_diagonal)

    # Pre-filled cells
    pre_filled = {
        (0, 1): 'g', (0, 2): 'f', (0, 4): 'e',
        (1, 1): 'f', (1, 2): 'a', (1, 4): 'c', (1, 6): 'b',
        (2, 3): 'c', (2, 4): 'd',
        (3, 2): 'c', (3, 3): 'd', (3, 4): 'b',
        (4, 0): 'e', (4, 1): 'c', (4, 6): 'a',
        (5, 1): 'd', (5, 3): 'g', (5, 4): 'f', (5, 5): 'a', (5, 6): 'e',
        (6, 2): 'g', (6, 3): 'f', (6, 5): 'e', (6, 6): 'c'
    }

    for (r, c), value in pre_filled.items():
        problem.addConstraint(lambda var, val=value: var == val, [(r, c)])

    # Solve the problem
    solution = problem.getSolution()

    if solution:
        # Format the solution
        grid = [['' for _ in range(7)] for _ in range(7)]
        for (r, c), letter in solution.items():
            grid[r][c] = letter
        return '\n'.join(','.join(row) for row in grid)
    else:
        return "No solution found"

# Solve the puzzle
solution = solve_puzzle()
print(f"<<<\n{solution}\n>>>")