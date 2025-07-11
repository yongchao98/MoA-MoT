# Read the input sum from the user
input_sum = input()

# The numbers are separated by '+', so we split the string by '+' to get a list of numbers
summands = input_sum.split('+')

# Sort the list of numbers in non-decreasing order
summands.sort()

# Join the sorted numbers back into a string with '+' as the separator
sorted_sum = '+'.join(summands)

# Print the final rearranged sum
print(sorted_sum)