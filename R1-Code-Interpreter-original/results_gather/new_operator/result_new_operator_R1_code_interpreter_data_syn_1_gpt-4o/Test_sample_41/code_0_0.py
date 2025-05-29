# Define the operation a@b
def operation(a, b):
    return 2 * a * b

# Calculate (2@10)@5
result = operation(operation(2, 10), 5)
print(result)