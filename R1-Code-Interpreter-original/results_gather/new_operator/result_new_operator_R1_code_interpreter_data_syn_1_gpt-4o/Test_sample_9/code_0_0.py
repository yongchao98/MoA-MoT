# Define the operation a & b
def operation(a, b):
    return (a + b) / b

# Calculate (1 & 3) & 3
result = operation(operation(1, 3), 3)
print(result)