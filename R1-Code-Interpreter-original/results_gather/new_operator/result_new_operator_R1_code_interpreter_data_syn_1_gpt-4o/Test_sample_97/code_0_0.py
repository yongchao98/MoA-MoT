# Define the operation a $ b
def operation(a, b):
    return (a + b) / b

# Calculate (3 $ 4) $ 10
result = operation(operation(3, 4), 10)
print(result)