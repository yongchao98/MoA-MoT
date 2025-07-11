# Step 1: Read the input sum as a string.
input_sum = input()

# Step 2: Split the string by the '+' character to get a list of numbers.
# For example, "3+2+1+4+2" becomes ['3', '2', '1', '4', '2']
summands = input_sum.split('+')

# Step 3: Sort the list of numbers. Python's default sort for strings is
# lexicographical, which works for single digits.
# ['3', '2', '1', '4', '2'] becomes ['1', '2', '2', '3', '4']
summands.sort()

# Step 4: Join the sorted numbers back into a string, separated by '+'.
# ['1', '2', '2', '3', '4'] becomes "1+2+2+3+4"
output_sum = "+".join(summands)

# Step 5: Print the final rearranged sum. The output format correctly shows each number in the equation.
print(output_sum)