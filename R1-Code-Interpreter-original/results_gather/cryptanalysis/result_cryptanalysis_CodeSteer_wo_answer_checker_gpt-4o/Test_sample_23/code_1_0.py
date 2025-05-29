from itertools import product

# Define possible values
numbers = '0123456789'
letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

# Initialize domains for each position
domains = {
    'pos1': set(numbers),
    'pos2': set(numbers),
    'pos3': set(letters),
    'pos4': set(letters)
}

# Apply constraints from each guess
def apply_constraints():
    # Guess: 98LZ
    domains['pos1'].intersection_update('9')
    domains['pos2'].difference_update('98')
    domains['pos3'].difference_update('LMNOPQRSTUVWXYZ')
    domains['pos4'].difference_update('LMNOPQRSTUVWXYZ')

    # Guess: 82EM
    domains['pos1'].difference_update('8')
    domains['pos2'].intersection_update('2')
    domains['pos3'].intersection_update('E')
    domains['pos4'].difference_update('M')

    # Guess: 36ZI
    domains['pos1'].difference_update('36')
    domains['pos2'].difference_update('36')
    domains['pos3'].difference_update('Z')
    domains['pos4'].difference_update('ZI')

    # Guess: 21HR
    domains['pos1'].intersection_update('2')
    domains['pos2'].difference_update('1')
    domains['pos3'].intersection_update('H')
    domains['pos4'].difference_update('R')

# Backtracking search
def backtrack(assignment={}):
    # If assignment is complete, return it
    if len(assignment) == 4:
        return assignment

    # Select an unassigned variable
    unassigned = [v for v in domains if v not in assignment]
    first = unassigned[0]

    # Try each value in the domain of the variable
    for value in domains[first]:
        # Assign value
        local_assignment = assignment.copy()
        local_assignment[first] = value

        # Check consistency
        if is_consistent(local_assignment):
            # Recurse
            result = backtrack(local_assignment)
            if result is not None:
                return result

    return None

# Check if the current assignment is consistent
def is_consistent(assignment):
    # Check if the assignment satisfies all constraints
    # In this case, we assume the constraints have been applied correctly
    return True

# Apply constraints
apply_constraints()

# Find the solution
solution = backtrack()

# Output the solution
if solution:
    password = [solution['pos1'], solution['pos2'], solution['pos3'], solution['pos4']]
    print(password)
else:
    print("No solution found")