# First, ensure you have the library installed:
# pip install python-constraint

from constraint import Problem, AllDifferentConstraint

def solve_puzzle():
    problem = Problem()

    # Define variables for each cell in the grid
    variables = [(r, c) for r in range(7) for c in range(7)]
    letters = 'abcdefg'
    problem.addVariables(variables, letters)

    # Add constraints for each row and column to have unique letters
    for i in range(7):
        problem.addConstraint(AllDifferentConstraint(), [(i, j) for j in range(7)])  # Row constraint
        problem.addConstraint(AllDifferentConstraint(), [(j, i) for j in range(7)])  # Column constraint

    # Add constraint for the minor diagonal to have the same letter
    minor_diagonal = [(i, 6 - i) for i in range(7)]
    for i in range(1, 7):
        problem.addConstraint(lambda a, b: a == b, (minor_diagonal[0], minor_diagonal[i]))

    # Add constraints for pre-filled cells
    pre_filled = {
        (0, 0): 'a', (0, 1): 'e', (0, 2): 'f', (0, 3): 'g', (0, 5): 'c', (0, 6): 'd',
        (1, 0): 'e', (1, 1): 'f', (1, 2): 'g', (1, 3): 'b', (1, 5): 'd', (1, 6): 'a',
        (2, 1): 'g', (2, 2): 'b', (2, 6): 'e',
        (3, 1): 'b', (3, 4): 'a', (3, 6): 'f',
        (4, 1): 'c', (4, 2): 'd', (4, 4): 'e', (4, 6): 'g',
        (5, 0): 'c', (5, 1): 'd', (5, 2): 'a', (5, 3): 'e', (5, 4): 'f', (5, 5): 'g', (5, 6): 'b',
        (6, 0): 'd', (6, 1): 'a', (6, 2): 'e', (6, 3): 'f', (6, 5): 'b'
    }
    for (r, c), letter in pre_filled.items():
        problem.addConstraint(lambda var, val=letter: var == val, [(r, c)])

    # Solve the problem
    solution = problem.getSolution()

    if solution:
        # Print the solution in the required format
        for r in range(7):
            print(','.join(solution[(r, c)] for c in range(7)))
    else:
        print("No solution found")

solve_puzzle()