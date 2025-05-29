# Define the operations
def operation_hash(a, b):
    return 4 * a * b

def operation_tilde(a, b):
    return 3 * a + 5 * b

def operation_percent(a, b):
    return 4 * a * b

# Calculate each step
step1 = operation_hash(3, 4)
step2 = operation_tilde(step1, 5)
step3 = operation_percent(step2, 6)

print(step3)