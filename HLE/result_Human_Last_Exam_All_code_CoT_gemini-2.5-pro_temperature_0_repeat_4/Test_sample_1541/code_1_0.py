# Read the input sum as a string from the user.
# For example, if the input is "3+2+1+4", the variable 's' will hold this string.
s = input()

# The summands are separated by '+'. We use the split() method to separate the string
# into a list of numbers. For "3+2+1+4", this creates the list ['3', '2', '1', '4'].
numbers = s.split('+')

# Sort the list of numbers. For single-digit numbers, lexicographical sorting
# (alphabetical) is the same as numerical sorting.
# The list ['3', '2', '1', '4'] becomes ['1', '2', '3', '4'].
numbers.sort()

# Join the elements of the sorted list back into a single string, with '+'
# as the separator between the numbers.
# The list ['1', '2', '3', '4'] becomes the string "1+2+3+4".
new_sum = '+'.join(numbers)

# Print the final rearranged sum.
print(new_sum)