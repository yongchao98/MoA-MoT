# Define the operations
def operation_hash(a, b):
    return 4 * a + 1 * b

def operation_at(a, b):
    return 2 * a * b

def operation_tilde(a, b):
    return 3 * a + 3 * b

# Calculate the expression (3#4)@5~6
result_hash = operation_hash(3, 4)
result_at = operation_at(result_hash, 5)
result_tilde = operation_tilde(result_at, 6)

print(result_hash, result_at, result_tilde)