# Define the operations
def operation_at(a, b):
    return (a + b) / b

def operation_percent(a, b):
    return 4 * a + 4 * b

def operation_hash(a, b):
    return 1 * a + 1 * b

# Calculate each step
result_at = operation_at(3, 4)
result_percent = operation_percent(result_at, 5)
result_hash = operation_hash(result_percent, 6)

# Print the final result
print(result_hash)