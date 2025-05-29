# Define the operations
def operation_hash(a, b):
    return (a + b) / b

def operation_tilde(a, b):
    return 4 * a + 5 * b

def operation_dollar(a, b):
    return (a + b) / b

# Calculate each part
result_hash = operation_hash(3, 4)
result_dollar = operation_dollar(5, 6)

# Calculate the final result
final_result = operation_tilde(result_hash, result_dollar)

print(final_result)