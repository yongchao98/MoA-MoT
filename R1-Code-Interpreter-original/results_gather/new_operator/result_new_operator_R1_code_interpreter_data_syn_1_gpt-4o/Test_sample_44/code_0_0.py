# Define the operations
def operation_dollar(a, b):
    return 1 * a * b

def operation_tilde(a, b):
    return 5 * a + 5 * b

def operation_ampersand(a, b):
    return 4 * a * b

# Calculate each step
result_dollar = operation_dollar(3, 4)
result_tilde = operation_tilde(result_dollar, 5)
result_ampersand = operation_ampersand(result_tilde, 6)

# Print the final result
print(result_ampersand)