# Read the input string from the user
s = input()

# Split the string by the '+' character to get a list of number strings
numbers = s.split('+')

# Sort the list of numbers. Lexicographical sorting works correctly for single digits.
numbers.sort()

# Join the sorted list of numbers back into a single string with '+' as a separator
new_sum = '+'.join(numbers)

# Print the final rearranged sum
print(new_sum)