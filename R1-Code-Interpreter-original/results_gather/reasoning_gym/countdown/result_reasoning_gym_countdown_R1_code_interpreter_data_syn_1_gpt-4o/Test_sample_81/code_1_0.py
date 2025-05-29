# Define the numbers
numbers = [23, 58, 5, 5]

# Function to evaluate expressions
def evaluate_expression():
    # Try different combinations of operations
    # Multiply 58 and 5, then add 23 and 5
    result1 = 58 * 5 + 23 + 5
    # Multiply 23 and 5, then multiply by 5 and add 58
    result2 = 23 * 5 * 5 + 58
    # Multiply 23 and 58, then subtract 5 and 5
    result3 = 23 * 58 - 5 - 5
    # Multiply 23 and 58, then add 5 and 5
    result4 = 23 * 58 + 5 + 5
    # Multiply 23 and 5, then add 58 and 5
    result5 = 23 * 5 + 58 + 5
    # Multiply 58 and 5, then add 23 and 5
    result6 = 58 * 5 + 23 + 5

    # Check which result equals 262
    if result1 == 262:
        return "58 * 5 + 23 + 5"
    elif result2 == 262:
        return "23 * 5 * 5 + 58"
    elif result3 == 262:
        return "23 * 58 - 5 - 5"
    elif result4 == 262:
        return "23 * 58 + 5 + 5"
    elif result5 == 262:
        return "23 * 5 + 58 + 5"
    elif result6 == 262:
        return "58 * 5 + 23 + 5"
    else:
        return "No solution found"

# Print the result
print(evaluate_expression())