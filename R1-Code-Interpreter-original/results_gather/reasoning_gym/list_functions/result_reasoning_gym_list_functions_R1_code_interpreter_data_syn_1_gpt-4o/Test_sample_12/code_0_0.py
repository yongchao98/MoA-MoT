# Define the input list
input_list = [76, 13, 35, 41]

# Initialize the output list with the first element unchanged
output_list = [input_list[0]]

# Apply the pattern to the rest of the elements
for i in range(1, len(input_list)):
    output_list.append(input_list[i] + i)

# Print the output list
print(output_list)