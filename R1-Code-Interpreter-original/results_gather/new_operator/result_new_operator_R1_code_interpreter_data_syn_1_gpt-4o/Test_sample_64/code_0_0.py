# Define the operations
def operation_hash(a, b):
    return 1 * a + 4 * b

def operation_at(a, b):
    return (a + b) / b

def operation_and(a, b):
    return (a + b) / b

# Calculate step by step
result_hash = operation_hash(3, 4)
result_at = operation_at(result_hash, 5)
result_and = operation_and(result_at, 6)

# Print the final result
print(result_and)