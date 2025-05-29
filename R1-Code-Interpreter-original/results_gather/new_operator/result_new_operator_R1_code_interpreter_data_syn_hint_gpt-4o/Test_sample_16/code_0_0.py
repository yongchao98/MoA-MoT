# Define the operations
def operation_at(a, b):
    return (a + b) / b

def operation_hash(a, b):
    return 2 * a + 3 * b

def operation_dollar(a, b):
    return (a + b) / b

# Calculate the expression (3@4)#5$6
result_at = operation_at(3, 4)
result_hash = operation_hash(result_at, 5)
result_dollar = operation_dollar(result_hash, 6)

# Print the final result
print(result_dollar)