# Define the operations
def operation_percent(a, b):
    return (a + b) / b

def operation_hash(a, b):
    return 5 * a + 1 * b

def operation_tilde(a, b):
    return 3 * a * b

# Calculate the expression (3%4)#5~6
result_percent = operation_percent(3, 4)
result_hash = operation_hash(result_percent, 5)
result_tilde = operation_tilde(result_hash, 6)

print(result_tilde)