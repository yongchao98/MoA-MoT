from constraint import Problem, AllDifferentConstraint

def solve_n_queens_with_constraints(pre_existing, blocked):
    problem = Problem()

    # Define variables for each row with possible column positions
    for row in range(8):
        if row in pre_existing:
            # If there's a pre-existing queen, fix its position
            problem.addVariable(row, [pre_existing[row]])
        else:
            # Otherwise, allow any column except blocked ones
            problem.addVariable(row, [col for col in range(8) if (row, col) not in blocked])

    # Add constraints for columns and diagonals
    problem.addConstraint(AllDifferentConstraint(), range(8))  # All columns must be different

    # Add diagonal constraints
    def diagonal_constraint(*args):
        for i in range(len(args)):
            for j in range(i + 1, len(args)):
                if abs(args[i] - args[j]) == abs(i - j):
                    return False
        return True

    problem.addConstraint(diagonal_constraint, range(8))

    # Solve the problem
    solutions = problem.getSolutions()

    if solutions:
        # Print the first solution found
        solution = solutions[0]
        positions = [f"{row} {solution[row]}" for row in range(8)]
        print("<<<" + ", ".join(positions) + ">>>")
    else:
        print("No solution exists")

def main():
    # Example pre-existing queens and blocked positions
    pre_existing = {}  # e.g., {0: 0} means a queen is already at (0, 0)
    blocked = set()    # e.g., {(1, 1)} means (1, 1) is blocked

    solve_n_queens_with_constraints(pre_existing, blocked)

main()