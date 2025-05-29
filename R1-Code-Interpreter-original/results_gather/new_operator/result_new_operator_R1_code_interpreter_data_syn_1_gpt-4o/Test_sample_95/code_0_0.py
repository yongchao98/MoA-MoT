# Define the operation a ~ b
def operation(a, b):
    return (a + b) / b

# Calculate (4 ~ 8) ~ 4
result = operation(operation(4, 8), 4)
print(result)