# Read the sum from the input
input_sum = input()

# Split the string by the '+' character to get a list of the numbers
summands = input_sum.split('+')

# Sort the list of numbers. Lexicographical sort works correctly for single digits.
summands.sort()

# Join the sorted list back into a string with '+' as the separator
output_sum = "+".join(summands)

# Print the final rearranged sum
print(output_sum)