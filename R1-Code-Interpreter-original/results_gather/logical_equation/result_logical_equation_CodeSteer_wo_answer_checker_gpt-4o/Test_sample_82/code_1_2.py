from itertools import permutations

# Given numbers
numbers = [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]

# Constraints as functions
def constraints(assignment):
    A, B, C, D, E, F, G, H, I, J, K, L, M = assignment
    return (
        C - E == -21 and
        J == 4 * C and
        H + K == 23 and
        D + L == 17 and
        C - A == -26 and
        G + M == 8 and
        F + G == 13 and
        H - G == 4 and
        A - C == 26 and
        I - J == -68 and
        F + M == 15
    )

# Function to solve the problem using backtracking
def solve():
    # Try all permutations of the numbers
    for perm in permutations(numbers):
        # Check if the current permutation satisfies all constraints
        if constraints(perm):
            return perm

# Find the solution
solution = solve()

# Print the solution
if solution:
    print(f"<<<{list(solution)}>>>")
else:
    print("No solution found.")