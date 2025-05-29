from itertools import permutations

# Possible values for the letters
values = [1, 2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50]

# Function to check if the current assignment satisfies all constraints
def satisfies_constraints(assignment):
    A, B, C, D, E, F, G, H, I, J, K, L, M = assignment
    
    return (A - G == -13 and
            B + C == 15 and
            H == 4.8 * B and
            I - B == -3 and
            L - J == -5 and
            B > A and
            I == 2.0 * F and
            C - H == -14 and
            B - I == 3 and
            A - F == 2)

# Backtracking function to find a valid assignment
def find_assignment(assignment, remaining_values):
    if len(assignment) == 13:
        if satisfies_constraints(assignment):
            return assignment
        return None
    
    for value in remaining_values:
        new_assignment = assignment + [value]
        new_remaining_values = remaining_values - {value}
        result = find_assignment(new_assignment, new_remaining_values)
        if result is not None:
            return result
    
    return None

# Start with an empty assignment and all values available
initial_assignment = []
remaining_values = set(values)

# Find a valid assignment
solution = find_assignment(initial_assignment, remaining_values)

# Print the solution if found
if solution:
    print(solution)