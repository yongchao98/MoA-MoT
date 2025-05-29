# Define the operation
def custom_operation(a, b):
    return 4 * a + 4 * b

# Calculate 3 $ 7
result_3_7 = custom_operation(3, 7)

# Calculate (3 $ 7) $ 1
final_result = custom_operation(result_3_7, 1)

# Print the final result
print(final_result)