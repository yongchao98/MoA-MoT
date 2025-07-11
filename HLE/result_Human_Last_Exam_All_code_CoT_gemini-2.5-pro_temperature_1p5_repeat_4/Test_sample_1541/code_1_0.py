# Read the input sum as a string
s = input()

# The summands are the numbers separated by '+'.
# We can use split('+') to get a list of these numbers as strings.
# For example, if s is "3+2+1+4+2", numbers will be ['3', '2', '1', '4', '2']
numbers = s.split('+')

# Sort the list of numbers. Since they are single-digit numbers,
# a standard lexicographical sort works correctly.
# For example, ['3', '2', '1', '4', '2'] becomes ['1', '2', '2', '3', '4']
numbers.sort()

# Join the sorted numbers back into a string, with '+' between them.
# For example, ['1', '2', '2', '3', '4'] becomes "1+2+2+3+4"
new_sum = "+".join(numbers)

# Print the final rearranged sum.
print(new_sum)