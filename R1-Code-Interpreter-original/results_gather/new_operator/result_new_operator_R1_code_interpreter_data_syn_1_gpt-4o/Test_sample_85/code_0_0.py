# Define the operation a & b
def operation(a, b):
    return (a + b) / b

# Calculate (7 & 5) & 4
result = operation(operation(7, 5), 4)
print(result)