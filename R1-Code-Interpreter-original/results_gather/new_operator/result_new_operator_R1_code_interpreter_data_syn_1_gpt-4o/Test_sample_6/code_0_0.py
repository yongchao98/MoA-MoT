# Define the operations
def operation_tilde(a, b):
    return 4 * a + 5 * b

def operation_hash(a, b):
    return (a + b) / b

def operation_ampersand(a, b):
    return 1 * a * b

# Calculate the expression (3~4)#5&6
result_tilde = operation_tilde(3, 4)
result_hash = operation_hash(result_tilde, 5)
result_ampersand = operation_ampersand(result_hash, 6)

print(result_ampersand)