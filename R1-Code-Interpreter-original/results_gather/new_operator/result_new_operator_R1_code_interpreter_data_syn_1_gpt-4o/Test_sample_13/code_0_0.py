# Define the operation a ~ b
def operation(a, b):
    return (a + b) / b

# Calculate (1 ~ 5) ~ 3
result = operation(operation(1, 5), 3)
print(result)