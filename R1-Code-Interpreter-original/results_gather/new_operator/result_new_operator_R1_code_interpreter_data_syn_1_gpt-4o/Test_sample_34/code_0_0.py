# Define the operations
def operation_percent(a, b):
    return 3 * a + 3 * b

def operation_dollar(a, b):
    return (a + b) / b

def operation_at(a, b):
    return 5 * a + 1 * b

# Calculate the expression (3%4)$5@6
result_percent = operation_percent(3, 4)
result_dollar = operation_dollar(result_percent, 5)
result_at = operation_at(result_dollar, 6)

# Print the final result
print(result_at)