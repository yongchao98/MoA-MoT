# Define the operation
def operation(a, b):
    return 3 * a * b

# Calculate (10 ~ 10) ~ 3
result = operation(operation(10, 10), 3)
print(result)