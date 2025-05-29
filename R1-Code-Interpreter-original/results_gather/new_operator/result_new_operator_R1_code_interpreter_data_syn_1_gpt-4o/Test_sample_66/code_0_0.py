# Define the operations
def at_operation(a, b):
    return 4 * a + 2 * b

def dollar_operation(a, b):
    return 5 * a + 1 * b

def percent_operation(a, b):
    return 4 * a + 4 * b

# Calculate the expression (3@4)$5%6
result = percent_operation(dollar_operation(at_operation(3, 4), 5), 6)
print(result)