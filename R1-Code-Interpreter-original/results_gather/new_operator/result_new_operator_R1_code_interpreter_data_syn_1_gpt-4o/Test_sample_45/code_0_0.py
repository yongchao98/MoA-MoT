# Define the operation a ~ b
def operation(a, b):
    return 3 * a * b

# Calculate 10 ~ 10
result1 = operation(10, 10)

# Calculate (10 ~ 10) ~ 3
result2 = operation(result1, 3)

# Print the final result
print(result2)