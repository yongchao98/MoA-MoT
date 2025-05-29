# Define the input list
input_list = [25, 100, 25, 1, 68]

# Apply the rule: increment the first and last elements by 1
output_list = input_list[:]
output_list[0] += 1
output_list[-1] += 1

# Print the output list
print(output_list)