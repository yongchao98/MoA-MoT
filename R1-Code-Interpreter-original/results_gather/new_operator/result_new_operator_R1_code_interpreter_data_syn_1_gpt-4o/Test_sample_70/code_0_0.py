# Define the operations
def operation_hash(a, b):
    return (a + b) / b

def operation_tilde(a, b):
    return 3 * a + 5 * b

def operation_ampersand(a, b):
    return 1 * a * b

# Calculate each part
result_hash = operation_hash(3, 4)
result_ampersand = operation_ampersand(5, 6)
final_result = operation_tilde(result_hash, 5)

# Print the final result
print(final_result)