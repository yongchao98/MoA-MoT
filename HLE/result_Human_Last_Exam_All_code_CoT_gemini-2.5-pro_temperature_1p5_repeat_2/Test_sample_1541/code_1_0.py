# Read the sum from the input as a string
s = input()

# Split the string by the '+' character to get a list of numbers (as strings)
numbers = s.split('+')

# Sort the list of numbers. For single-digit numbers, default string sort is sufficient.
numbers.sort()

# Join the sorted numbers back into a string with '+' as the separator
sorted_sum = "+".join(numbers)

# Print the final rearranged sum
print(sorted_sum)