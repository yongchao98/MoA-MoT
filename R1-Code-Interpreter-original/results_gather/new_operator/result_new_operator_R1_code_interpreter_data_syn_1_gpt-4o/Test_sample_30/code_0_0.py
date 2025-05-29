# Define the operations
def operation_tilde(a, b):
    return 3 * a + 1 * b

def operation_and(a, b):
    return 3 * a * b

def operation_hash(a, b):
    return 1 * a * b

# Calculate the expression (3~4)&5#6
result_tilde = operation_tilde(3, 4)
result_and = operation_and(result_tilde, 5)
result_hash = operation_hash(result_and, 6)

# Print the final result
print(result_hash)