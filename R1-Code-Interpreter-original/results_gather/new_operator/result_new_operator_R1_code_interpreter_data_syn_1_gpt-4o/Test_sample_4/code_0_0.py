# Define the operations
def operation_percent(a, b):
    return 4 * a + 2 * b

def operation_hash(a, b):
    return (a + b) / b

def operation_ampersand(a, b):
    return (a + b) / b

# Calculate the expression (3%4)#5&6
result_percent = operation_percent(3, 4)
result_hash = operation_hash(result_percent, 5)
result_ampersand = operation_ampersand(result_hash, 6)

print(result_percent, result_hash, result_ampersand)