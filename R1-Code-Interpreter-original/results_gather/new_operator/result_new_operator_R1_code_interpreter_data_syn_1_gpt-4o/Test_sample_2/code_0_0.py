# Define the operations
def operation_and(a, b):
    return 1 * a * b

def operation_tilde(a, b):
    return 3 * a * b

def operation_at(a, b):
    return 2 * a + 2 * b

# Calculate the expression (3&4)~5@6
result = operation_at(operation_tilde(operation_and(3, 4), 5), 6)
print(result)