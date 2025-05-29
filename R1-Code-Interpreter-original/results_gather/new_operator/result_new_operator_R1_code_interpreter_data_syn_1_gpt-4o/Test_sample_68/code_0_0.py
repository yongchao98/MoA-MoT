# Define the operations
def operation_at(a, b):
    return 4 * a + 2 * b

def operation_and(a, b):
    return (a + b) / b

def operation_hash(a, b):
    return (a + b) / b

# Calculate the expression (3@4)&5#6
result = operation_hash(operation_and(operation_at(3, 4), 5), 6)
print(result)