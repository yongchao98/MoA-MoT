# Define the operation a $ b
def operation(a, b):
    return (a + b) / b

# Calculate (10 $ 3) $ 8
result = operation(operation(10, 3), 8)
print(result)