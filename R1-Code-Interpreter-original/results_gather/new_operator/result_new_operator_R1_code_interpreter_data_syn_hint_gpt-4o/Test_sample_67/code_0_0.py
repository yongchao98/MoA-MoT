# Define the operation
def operation(a, b):
    return 5 * a * b

# Calculate (1 ~ 10) ~ 9
result = operation(operation(1, 10), 9)
print(result)