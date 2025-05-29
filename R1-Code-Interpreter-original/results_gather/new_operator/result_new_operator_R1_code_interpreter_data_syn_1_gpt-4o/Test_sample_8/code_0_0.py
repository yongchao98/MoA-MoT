# Define the operations
def operation_percent(a, b):
    return 3 * a * b

def operation_dollar(a, b):
    return 5 * a + 4 * b

def operation_hash(a, b):
    return 2 * a + 2 * b

# Calculate the expression (3%4)$5#6
result_percent = operation_percent(3, 4)
result_dollar = operation_dollar(result_percent, 5)
result_hash = operation_hash(result_dollar, 6)

print(result_hash)