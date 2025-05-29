# Define the operations
def operation_hash(a, b):
    return 1 * a + 5 * b

def operation_tilde(a, b):
    return 5 * a * b

def operation_dollar(a, b):
    return 1 * a * b

# Calculate each part
result_hash = operation_hash(3, 4)
result_dollar = operation_dollar(5, 6)
final_result = operation_tilde(result_hash, result_dollar)

# Print the final result
print(final_result)