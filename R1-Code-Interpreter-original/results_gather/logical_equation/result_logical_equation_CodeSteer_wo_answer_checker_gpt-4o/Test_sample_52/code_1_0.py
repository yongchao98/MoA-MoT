from itertools import permutations

# Define the possible values for each letter
possible_values = [3, 9, 16, 20, 39, 48, 75, 80, 121, 150, 225]

# Function to check if a given assignment satisfies all constraints
def satisfies_constraints(assignment):
    A, B, C, D, E, F, G, H, I, J, K = assignment
    
    # Check all constraints
    return (
        D - A == 66 and
        F - G == -145 and
        B - G == -209 and
        E == 3.0 * B and
        I - A == 112 and
        A + H == 159 and
        D + H == 225 and
        E == 2.4 * J and
        A - H == -141
    )

# Backtracking function to find a valid solution
def find_solution():
    for perm in permutations(possible_values):
        if satisfies_constraints(perm):
            return perm
    return None

# Find and print the solution
solution = find_solution()
if solution:
    print(f"<<<{list(solution)}>>>")
else:
    print("No valid solution found that matches the possible values.")