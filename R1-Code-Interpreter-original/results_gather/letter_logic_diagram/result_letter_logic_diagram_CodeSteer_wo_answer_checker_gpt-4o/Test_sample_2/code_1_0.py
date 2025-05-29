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
    minor_diagonal = [(0, 6), (1, 5), (2, 4), (3, 3), (4, 2), (5, 1), (6, 0)]
    problem.addConstraint(lambda *args: len(set(args)) == 1, minor_diagonal)

    # Pre-filled cells
    pre_filled = {
        (0, 1): 'b', (0, 3): 'f', (0, 4): 'g', (0, 6): 'c',
        (1, 0): 'b', (1, 1): 'e', (1, 4): 'a', (1, 5): 'c',
        (2, 1): 'f', (2, 2): 'g', (2, 3): 'a', (2, 4): 'c', (2, 6): 'b',
        (3, 0): 'f', (3, 2): 'a', (3, 4): 'd', (3, 6): 'e',
        (4, 0): 'g', (4, 3): 'd', (4, 5): 'e', (4, 6): 'f',
        (5, 0): 'a', (5, 5): 'f', (5, 6): 'g',
        (6, 0): 'c', (6, 1): 'd', (6, 3): 'e', (6, 5): 'g', (6, 6): 'a'
    }

    for (r, c), value in pre_filled.items():
        problem.addConstraint(lambda var, val=value: var == val, [(r, c)])

    # Solve the problem
    solution = problem.getSolution()

    if solution:
        # Format the output
        result = []
        for r in range(7):
            row = [solution[(r, c)] for c in range(7)]
            result.append(','.join(row))
        return f"<<<\n" + '\n'.join(result) + "\n>>>"
    else:
        return "No solution found"

# Solve the puzzle
print(solve_puzzle_with_csp())