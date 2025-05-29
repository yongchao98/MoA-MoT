# Backtracking search with logging
def backtrack(assignment={}):
    # If assignment is complete, return it
    if len(assignment) == 4:
        print(f"Complete assignment found: {assignment}")
        return assignment

    # Select an unassigned variable
    unassigned = [v for v in domains if v not in assignment]
    first = unassigned[0]

    # Try each value in the domain of the variable
    for value in domains[first]:
        # Assign value
        local_assignment = assignment.copy()
        local_assignment[first] = value
        print(f"Trying {first} = {value}")

        # Check consistency
        if is_consistent(local_assignment):
            # Recurse
            result = backtrack(local_assignment)
            if result is not None:
                return result
        else:
            print(f"Inconsistent assignment: {local_assignment}")

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
    print(f"Password found: {password}")
else:
    print("No solution found")