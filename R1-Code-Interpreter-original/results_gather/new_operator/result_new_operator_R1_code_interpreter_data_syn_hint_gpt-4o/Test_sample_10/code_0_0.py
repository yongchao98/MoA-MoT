# Define the operations
def operation_tilde(a, b):
    return (a + b) / b

def operation_hash(a, b):
    return 5 * a + 5 * b

def operation_at(a, b):
    return 1 * a * b

# Calculate the expression step by step
result_tilde = operation_tilde(3, 4)
result_hash = operation_hash(result_tilde, 5)
result_at = operation_at(result_hash, 6)

# Print the final result
print(result_at)