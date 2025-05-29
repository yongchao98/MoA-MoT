# Define the operation a ~ b
def operation(a, b):
    return (a + b) / b

# Calculate (6 ~ 1) ~ 9
result = operation(operation(6, 1), 9)
print(result)