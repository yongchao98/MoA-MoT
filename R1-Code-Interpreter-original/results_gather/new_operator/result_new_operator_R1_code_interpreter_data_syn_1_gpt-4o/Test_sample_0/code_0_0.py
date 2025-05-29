# Define the operations
def operation_tilde(a, b):
    return 1 * a * b

def operation_and(a, b):
    return 5 * a + 5 * b

def operation_at(a, b):
    return 3 * a * b

# Calculate each step
result_tilde = operation_tilde(3, 4)
result_and = operation_and(result_tilde, 5)
result_at = operation_at(5, 6)

# Final calculation
final_result = operation_at(result_and, result_at)

# Print the final result
print(final_result)