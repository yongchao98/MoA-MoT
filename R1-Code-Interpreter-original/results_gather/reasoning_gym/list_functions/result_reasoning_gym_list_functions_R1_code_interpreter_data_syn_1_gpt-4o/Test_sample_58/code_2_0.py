# Calculate the sum of the input list
input_list = [4, 16, 65, 37]
sum_of_elements = sum(input_list)

# Calculate the sum of the digits of the sum
sum_of_digits = sum(int(digit) for digit in str(sum_of_elements))

# Print the sum of the digits
print(sum_of_digits)