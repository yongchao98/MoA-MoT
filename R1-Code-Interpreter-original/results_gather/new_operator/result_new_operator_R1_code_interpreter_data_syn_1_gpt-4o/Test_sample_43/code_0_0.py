# Define the custom operation
def custom_mod(a, b):
    return 2 * a + 2 * b

# Step 1: Calculate 6%9
result_1 = custom_mod(6, 9)

# Step 2: Calculate result_1%5
final_result = custom_mod(result_1, 5)

# Print the final result
print(final_result)