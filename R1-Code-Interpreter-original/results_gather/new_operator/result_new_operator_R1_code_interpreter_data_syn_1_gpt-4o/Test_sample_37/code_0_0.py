# Define the operation a ~ b
def operation(a, b):
    return 3 * a + 2 * b

# Calculate 10 ~ 4
result1 = operation(10, 4)

# Calculate (10 ~ 4) ~ 3
result2 = operation(result1, 3)

# Print the final result
print(result2)