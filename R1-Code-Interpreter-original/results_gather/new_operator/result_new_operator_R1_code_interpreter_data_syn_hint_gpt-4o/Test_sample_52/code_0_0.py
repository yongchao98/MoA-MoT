# Define the operations
def operation_at(a, b):
    return 5 * a + 2 * b

def operation_percent(a, b):
    return 4 * a + 1 * b

def operation_and(a, b):
    return (a + b) / b

# Calculate the expression (3@4)%5&6
result_at = operation_at(3, 4)
result_percent = operation_percent(result_at, 5)
result_and = operation_and(result_percent, 6)

print(result_and)