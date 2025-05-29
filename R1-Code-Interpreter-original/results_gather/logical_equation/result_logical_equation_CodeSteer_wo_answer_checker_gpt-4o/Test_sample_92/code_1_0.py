def is_valid_partial(assignment, index):
    A, B, C, D, E, F, G, H, I, J, K, L, M = assignment
    constraints = [
        (A == 3.0 * L, index >= 0),
        (A + M == 39, index >= 0),
        (L + M == 29, index >= 0),
        (B == 3.2 * L, index >= 1),
        (A == 1.5 * I, index >= 0),
        (C - B == 80, index >= 2),
        (E + K == 43, index >= 4),
        (C > D, index >= 2),
        (E == 2.4 * A, index >= 4),
        (C + D == 124, index >= 2),
        (H - A == -12, index >= 7),
        (F + G == 47, index >= 5)
    ]
    return all(cond for cond, active in constraints if active)

def forward_checking(assignment, remaining_values, index):
    if index == 13:
        return assignment if is_valid_partial(assignment, index) else None

    for value in sorted(remaining_values):
        new_assignment = assignment + [value]
        if is_valid_partial(new_assignment, index):
            new_remaining_values = remaining_values - {value}
            result = forward_checking(new_assignment, new_remaining_values, index + 1)
            if result is not None:
                return result

    return None

# Define the possible values for each letter
possible_values = {2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96}

# Start the forward checking process
solution = forward_checking([], possible_values, 0)

# Print the solution in the required format
if solution:
    print(f"<<<{solution}>>>")
else:
    print("No solution found.")