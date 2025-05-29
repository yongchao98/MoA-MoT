# Define the operation a#b
def operation(a, b):
    return (a + b) / b

# Calculate (3#2)#6
result = operation(operation(3, 2), 6)
print(result)