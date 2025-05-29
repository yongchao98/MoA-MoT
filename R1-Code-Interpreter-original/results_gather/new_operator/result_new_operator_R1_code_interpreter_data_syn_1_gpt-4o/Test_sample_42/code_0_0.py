# Define the operations
def operation_and(a, b):
    return 4 * a + 3 * b

def operation_hash(a, b):
    return 3 * a + 2 * b

def operation_tilde(a, b):
    return 1 * a + 4 * b

# Calculate the expression (3&4)#5~6
result_and = operation_and(3, 4)
result_hash = operation_hash(result_and, 5)
result_tilde = operation_tilde(result_hash, 6)

# Print the final result
print(result_tilde)