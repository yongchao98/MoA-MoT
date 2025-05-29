# Define the operation a#b
def operation(a, b):
    return 3*a + 5*b

# Calculate 3#2
result_3_2 = operation(3, 2)

# Calculate (3#2)#5
final_result = operation(result_3_2, 5)

# Print the final result
print(final_result)