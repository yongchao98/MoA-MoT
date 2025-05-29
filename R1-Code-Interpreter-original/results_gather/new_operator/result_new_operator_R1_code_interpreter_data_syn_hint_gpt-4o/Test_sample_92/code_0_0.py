# Define the operations
def operation_tilde(a, b):
    return 5 * a + 4 * b

def operation_percent(a, b):
    return 2 * a * b

def operation_at(a, b):
    return 4 * a + 4 * b

# Calculate the expression (3~4)%5@6
result_tilde = operation_tilde(3, 4)
result_percent = operation_percent(result_tilde, 5)
result_at = operation_at(result_percent, 6)

# Print the final result
print(result_at)