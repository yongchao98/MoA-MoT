from itertools import permutations

# Define the possible values for each variable
possible_values = [3, 9, 16, 20, 39, 48, 75, 80, 121, 150, 225]

# Define the constraints as functions
def constraints(assignment):
    C, K, J, D, G, E, H, F, I, B = assignment['C'], assignment['K'], assignment['J'], assignment['D'], assignment['G'], assignment['E'], assignment['H'], assignment['F'], assignment['I'], assignment['B']
    return (
        C - K == -17 and
        J == 2 * D and
        G - K == 205 and
        E + H == 119 and
        C - J == -147 and
        F + K == 36 and
        F + G == 241 and
        H - F == 64 and
        I - B == -39 and
        D + K == 95 and
        B + J == 198
    )

# Backtracking function
def backtrack(assignment, variables):
    if len(assignment) == len(variables):
        return assignment

    # Select the next variable to assign
    unassigned_vars = [v for v in variables if v not in assignment]
    var = unassigned_vars[0]

    # Try each possible value for the variable
    for value in possible_values:
        if value not in assignment.values():
            assignment[var] = value
            if constraints(assignment):
                result = backtrack(assignment, variables)
                if result:
                    return result
            del assignment[var]

    return None

# Define the variables
variables = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K']

# Start the backtracking algorithm
solution = backtrack({}, variables)

# Print the result in the required format
if solution:
    result = [solution[var] for var in sorted(variables)]
    print(f"<<<{result}>>>")