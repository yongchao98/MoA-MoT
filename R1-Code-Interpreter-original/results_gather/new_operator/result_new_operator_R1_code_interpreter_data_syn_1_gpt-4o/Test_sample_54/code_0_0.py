# Define the operations
def operation_tilde(a, b):
    return (a + b) / b

def operation_percent(a, b):
    return (a + b) / b

def operation_dollar(a, b):
    return 2 * a + 4 * b

# Calculate step by step
result_tilde = operation_tilde(3, 4)
result_percent = operation_percent(result_tilde, 5)
result_dollar = operation_dollar(result_percent, 6)

# Print the final result
print(result_dollar)