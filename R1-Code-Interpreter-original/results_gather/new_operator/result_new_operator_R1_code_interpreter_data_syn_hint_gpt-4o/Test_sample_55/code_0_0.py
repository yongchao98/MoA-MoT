# Define the custom operation
def custom_mod(a, b):
    return 4 * a + 4 * b

# Step 1: Calculate 2 % 8
result_1 = custom_mod(2, 8)

# Step 2: Calculate (result_1) % 5
final_result = custom_mod(result_1, 5)

print(final_result)