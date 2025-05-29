# Define the operation a & b
def custom_and(a, b):
    return 4 * a + 2 * b

# Calculate 5 & 3
result_5_and_3 = custom_and(5, 3)

# Calculate (5 & 3) & 6
final_result = custom_and(result_5_and_3, 6)

# Print the final result
print(final_result)