from itertools import permutations

# Define the possible values for each variable
possible_values = [3, 9, 16, 20, 39, 48, 75, 80, 121, 150, 225]

# Define the constraints as functions
def check_constraints(assignment, var):
    # Check constraints involving the newly assigned variable
    if var == 'C':
        if 'K' in assignment and assignment['C'] - assignment['K'] != -17:
            return False
        if 'J' in assignment and assignment['C'] - assignment['J'] != -147:
            return False
    if var == 'K':
        if 'C' in assignment and assignment['C'] - assignment['K'] != -17:
            return False
        if 'G' in assignment and assignment['G'] - assignment['K'] != 205:
            return False
        if 'F' in assignment and assignment['F'] + assignment['K'] != 36:
            return False
        if 'D' in assignment and assignment['D'] + assignment['K'] != 95:
            return False
    if var == 'J':
        if 'D' in assignment and assignment['J'] != 2 * assignment['D']:
            return False
        if 'C' in assignment and assignment['C'] - assignment['J'] != -147:
            return False
        if 'B' in assignment and assignment['B'] + assignment['J'] != 198:
            return False
    if var == 'D':
        if 'J' in assignment and assignment['J'] != 2 * assignment['D']:
            return False
        if 'K' in assignment and assignment['D'] + assignment['K'] != 95:
            return False
    if var == 'G':
        if 'K' in assignment and assignment['G'] - assignment['K'] != 205:
            return False
        if 'F' in assignment and assignment['F'] + assignment['G'] != 241:
            return False
    if var == 'E':
        if 'H' in assignment and assignment['E'] + assignment['H'] != 119:
            return False
    if var == 'H':
        if 'E' in assignment and assignment['E'] + assignment['H'] != 119:
            return False
        if 'F' in assignment and assignment['H'] - assignment['F'] != 64:
            return False
    if var == 'F':
        if 'K' in assignment and assignment['F'] + assignment['K'] != 36:
            return False
        if 'G' in assignment and assignment['F'] + assignment['G'] != 241:
            return False
        if 'H' in assignment and assignment['H'] - assignment['F'] != 64:
            return False
    if var == 'I':
        if 'B' in assignment and assignment['I'] - assignment['B'] != -39:
            return False
    if var == 'B':
        if 'I' in assignment and assignment['I'] - assignment['B'] != -39:
            return False
        if 'J' in assignment and assignment['B'] + assignment['J'] != 198:
            return False
    return True

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
            if check_constraints(assignment, var):
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