def operation(a, b):
    return 2 * a + 4 * b

# Calculate 1@4
result1 = operation(1, 4)

# Calculate (1@4)@2
final_result = operation(result1, 2)

print(final_result)