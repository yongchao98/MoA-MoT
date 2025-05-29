# Define the operations
def operation_hash(a, b):
    return (a + b) / b

def operation_dollar(a, b):
    return 2 * a + 5 * b

def operation_percent(a, b):
    return (a + b) / b

# Calculate step by step
result_hash = operation_hash(3, 4)
result_dollar = operation_dollar(result_hash, 5)
result_percent = operation_percent(result_dollar, 6)

# Print the final result
print(result_percent)