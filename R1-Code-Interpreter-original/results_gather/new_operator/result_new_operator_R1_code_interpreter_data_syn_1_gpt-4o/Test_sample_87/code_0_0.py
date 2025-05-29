# Define the operation a & b
def operation(a, b):
    return (a + b) / b

# Calculate (3 & 5) & 2
result = operation(operation(3, 5), 2)
print(result)