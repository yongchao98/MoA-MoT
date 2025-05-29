# Define the operations
def operation_hash(a, b):
    return 2 * a * b

def operation_and(a, b):
    return 3 * a * b

def operation_at(a, b):
    return 4 * a * b

# Calculate the expression (3#4)&5@6
result_hash = operation_hash(3, 4)
result_and = operation_and(result_hash, 5)
result_at = operation_at(result_and, 6)

print(result_at)