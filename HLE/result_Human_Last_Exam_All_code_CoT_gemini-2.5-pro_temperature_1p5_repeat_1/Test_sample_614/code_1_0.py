# The riddle describes a specific sequence of five natural numbers.
# Based on a step-by-step logical deduction, the final sequence has been determined.
# This script will print each number of the solved sequence, which represents the solution
# to the "equation" proposed by the riddle.

solution_sequence = [1, 2, 3, 5, 4]

# The puzzle requests the numbers in the correct order.
# The following line prints each number from the solution list, separated by a comma and space.
print(*solution_sequence, sep=", ")