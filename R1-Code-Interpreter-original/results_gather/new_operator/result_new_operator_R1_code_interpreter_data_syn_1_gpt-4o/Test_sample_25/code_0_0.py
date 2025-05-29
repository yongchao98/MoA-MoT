# Define the operation a@b
def at_operation(a, b):
    return 1 * a * b

# Calculate (6@9)@4
result = at_operation(at_operation(6, 9), 4)
print(result)