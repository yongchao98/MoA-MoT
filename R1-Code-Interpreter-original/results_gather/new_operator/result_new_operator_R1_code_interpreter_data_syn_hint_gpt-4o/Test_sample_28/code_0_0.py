# Define the operations
def operation_tilde(a, b):
    return 1 * a * b

def operation_dollar(a, b):
    return (a + b) / b

def operation_at(a, b):
    return (a + b) / b

# Calculate the expression step by step
result_tilde = operation_tilde(3, 4)
result_dollar = operation_dollar(result_tilde, 5)
result_at = operation_at(result_dollar, 6)

# Print the final result
print(result_at)