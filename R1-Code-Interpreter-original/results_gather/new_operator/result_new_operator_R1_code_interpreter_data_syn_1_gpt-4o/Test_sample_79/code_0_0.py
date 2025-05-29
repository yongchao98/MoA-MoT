# Define the operation a ~ b
def operation(a, b):
    return 2 * a * b

# Calculate (1 ~ 9) ~ 1
result = operation(operation(1, 9), 1)
print(result)