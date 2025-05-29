def is_valid(assignment):
    """Check if the current assignment satisfies all constraints."""
    A, B, C, D, E, F, G, H, I, J, K, L, M = assignment

    return (
        E + F == 8 and
        C - D == 14 and
        A + M == 35 and
        D + K == 3 and
        G + J == 60 and
        H + I == 60 and
        A - K == 5 and
        G - H == 26 and
        G - F == 45 and
        B > F
    )

def backtrack(assignment, remaining_values):
    """Backtracking algorithm to find a valid assignment."""
    if len(assignment) == 13:
        if is_valid(assignment):
            return assignment
        return None

    for value in remaining_values:
        new_assignment = assignment + [value]
        new_remaining_values = remaining_values - {value}

        # Constraint propagation: Check if partial assignment is promising
        if len(new_assignment) >= 2:
            if not is_valid_partial(new_assignment):
                continue

        result = backtrack(new_assignment, new_remaining_values)
        if result is not None:
            return result

    return None

def is_valid_partial(assignment):
    """Check if the partial assignment is promising."""
    # Partial checks based on the constraints
    if len(assignment) >= 2:
        E, F = assignment[:2]
        if E + F != 8:
            return False
    if len(assignment) >= 4:
        C, D = assignment[2:4]
        if C - D != 14:
            return False
    if len(assignment) >= 6:
        A, M = assignment[4:6]
        if A + M != 35:
            return False
    if len(assignment) >= 8:
        D, K = assignment[6:8]
        if D + K != 3:
            return False
    if len(assignment) >= 10:
        G, J = assignment[8:10]
        if G + J != 60:
            return False
    if len(assignment) >= 12:
        H, I = assignment[10:12]
        if H + I != 60:
            return False
    if len(assignment) >= 13:
        A, K = assignment[4], assignment[7]
        if A - K != 5:
            return False
    return True

# Define the possible values
possible_values = {1, 2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50}

# Start the backtracking process
solution = backtrack([], possible_values)

# Print the solution if found
if solution:
    print(f"<<<{solution}>>>")
else:
    print("No valid solution found.")