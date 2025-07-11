# Read the input sum as a string
s = input()

# Split the string by the '+' character to get a list of numbers
numbers = s.split('+')

# Sort the list of numbers. Since they are single-digit strings,
# default string sorting works correctly.
numbers.sort()

# Join the sorted numbers back into a string with '+' as the separator
sorted_sum = '+'.join(numbers)

# Print the final rearranged sum
print(sorted_sum)