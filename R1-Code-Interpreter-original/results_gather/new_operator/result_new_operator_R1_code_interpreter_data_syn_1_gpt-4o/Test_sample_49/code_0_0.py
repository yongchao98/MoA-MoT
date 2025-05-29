# Define the operation
def dollar_operation(a, b):
    return (a + b) / b

# Calculate (3$10)$7
result = dollar_operation(dollar_operation(3, 10), 7)
print(result)