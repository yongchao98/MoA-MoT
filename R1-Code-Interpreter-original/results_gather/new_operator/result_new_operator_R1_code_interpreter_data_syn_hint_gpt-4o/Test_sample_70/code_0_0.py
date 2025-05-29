# Define the operations
def operation_hash(a, b):
    return (a + b) / b

def operation_tilde(a, b):
    return 3 * a + 5 * b

def operation_and(a, b):
    return 1 * a * b

# Calculate each part
result_hash = operation_hash(3, 4)
result_and = operation_and(5, 6)

# Calculate the final result
final_result = operation_tilde(result_hash, result_and)

print(final_result)