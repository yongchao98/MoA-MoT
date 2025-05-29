# Define the operations
def operation_hash(a, b):
    return 2 * a + 1 * b

def operation_tilde(a, b):
    return 1 * a + 3 * b

def operation_percent(a, b):
    return 5 * a * b

# Calculate the expression (3#4)~5%6
result_hash = operation_hash(3, 4)
result_tilde = operation_tilde(result_hash, 5)
result_percent = operation_percent(result_tilde, 6)

print(result_percent)