# Trying different combinations to reach 977
def find_expression():
    numbers = [11, 50, 92, 85]
    target = 977

    # Try multiplying 92 and 85
    result1 = 92 * 85
    # Check if adding or subtracting the remaining numbers gives the target
    if result1 + 11 + 50 == target:
        return f"92 * 85 + 11 + 50"
    if result1 + 11 - 50 == target:
        return f"92 * 85 + 11 - 50"
    if result1 - 11 + 50 == target:
        return f"92 * 85 - 11 + 50"
    if result1 - 11 - 50 == target:
        return f"92 * 85 - 11 - 50"

    # Try other combinations if needed
    return "No solution found with current approach"

print(find_expression())