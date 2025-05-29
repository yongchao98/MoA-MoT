# Trying different combinations to reach 977
def find_expression():
    numbers = [11, 50, 92, 85]
    target = 977

    # Try multiplying 50 and 11
    result1 = 50 * 11
    if result1 + 92 + 85 == target:
        return f"50 * 11 + 92 + 85"
    if result1 + 92 - 85 == target:
        return f"50 * 11 + 92 - 85"
    if result1 - 92 + 85 == target:
        return f"50 * 11 - 92 + 85"
    if result1 - 92 - 85 == target:
        return f"50 * 11 - 92 - 85"

    # Try multiplying 92 and 11
    result2 = 92 * 11
    if result2 + 50 + 85 == target:
        return f"92 * 11 + 50 + 85"
    if result2 + 50 - 85 == target:
        return f"92 * 11 + 50 - 85"
    if result2 - 50 + 85 == target:
        return f"92 * 11 - 50 + 85"
    if result2 - 50 - 85 == target:
        return f"92 * 11 - 50 - 85"

    # Try multiplying 85 and 11
    result3 = 85 * 11
    if result3 + 50 + 92 == target:
        return f"85 * 11 + 50 + 92"
    if result3 + 50 - 92 == target:
        return f"85 * 11 + 50 - 92"
    if result3 - 50 + 92 == target:
        return f"85 * 11 - 50 + 92"
    if result3 - 50 - 92 == target:
        return f"85 * 11 - 50 - 92"

    return "No solution found with current approach"

print(find_expression())