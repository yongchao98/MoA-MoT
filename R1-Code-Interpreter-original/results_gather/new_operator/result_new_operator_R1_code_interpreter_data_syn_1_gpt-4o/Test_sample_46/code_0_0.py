# Define the operations
def operation_and(a, b):
    return 3 * a * b

def operation_at(a, b):
    return (a + b) / b

def operation_percent(a, b):
    return (a + b) / b

# Calculate the expression (3&4)@5%6
result_and = operation_and(3, 4)
result_at = operation_at(result_and, 5)
result_percent = operation_percent(result_at, 6)

# Print the final result
print(result_percent)