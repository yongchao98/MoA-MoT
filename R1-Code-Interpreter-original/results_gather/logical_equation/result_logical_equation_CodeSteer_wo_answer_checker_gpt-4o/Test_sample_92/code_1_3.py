def is_valid(assignment):
    A, B, C, D, E, F, G, H, I, J, K, L, M = assignment
    return (
        A == 3.0 * L and
        A + M == 39 and
        L + M == 29 and
        B == 3.2 * L and
        A == 1.5 * I and
        C - B == 80 and
        E + K == 43 and
        C > D and
        E == 2.4 * A and
        C + D == 124 and
        H - A == -12 and
        F + G == 47
    )

def backtrack(assignment, remaining_values):
    if len(assignment) == 13:
        if is_valid(assignment):
            return assignment
        return None

    for value in remaining_values:
        new_assignment = assignment + [value]
        new_remaining_values = remaining_values - {value}
        result = backtrack(new_assignment, new_remaining_values)
        if result is not None:
            return result

    return None

# Define the possible values for each letter
possible_values = {2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96}

# Start the backtracking process
solution = backtrack([], possible_values)

# Print the solution in the required format
if solution:
    print(f"<<<{solution}>>>")
else:
    print("No solution found.")