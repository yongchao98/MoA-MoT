from constraint import Problem, AllDifferentConstraint

def solve_puzzle_with_csp():
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
    minor_diagonal = [(r, 6 - r) for r in range(7)]
    problem.addConstraint(lambda *args: len(set(args)) == 1, minor_diagonal)

    # Pre-filled cells
    pre_filled = {
        (0, 4): 'g',
        (2, 0): 'e', (2, 2): 'g', (2, 3): 'c', (2, 6): 'b',
        (3, 1): 'g', (3, 3): 'd', (3, 4): 'f', (3, 6): 'e',
        (4, 0): 'g', (4, 1): 'c', (4, 4): 'b', (4, 6): 'a',
        (5, 0): 'c', (5, 1): 'd', (5, 4): 'e', (5, 6): 'g',
        (6, 1): 'f', (6, 3): 'e'
    }
    for (r, c), value in pre_filled.items():
        problem.addConstraint(lambda var, val=value: var == val, [(r, c)])

    # Solve the problem
    solution = problem.getSolution()

    # Format the output
    if solution:
        grid = [['' for _ in range(7)] for _ in range(7)]
        for (r, c), value in solution.items():
            grid[r][c] = value
        return '\n'.join(','.join(row) for row in grid)
    else:
        return "No solution found"

# Solve the puzzle
solution = solve_puzzle_with_csp()
print(f"<<<\n{solution}\n>>>")