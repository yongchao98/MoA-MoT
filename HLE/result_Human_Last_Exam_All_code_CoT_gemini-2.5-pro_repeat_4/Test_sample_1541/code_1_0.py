# Read the input sum as a string
s = input()

# Split the string by '+' to get a list of numbers (as strings).
# For example, "3+2+1+4+1" becomes ['3', '2', '1', '4', '1'].
numbers = s.split('+')

# Sort the list of strings. Since all numbers are single digits,
# lexicographical (alphabetical) sorting is equivalent to numerical sorting.
# ['3', '2', '1', '4', '1'] becomes ['1', '1', '2', '3', '4'].
numbers.sort()

# Join the sorted list of strings back into a single string,
# separated by '+'.
# ['1', '1', '2', '3', '4'] becomes "1+1+2+3+4".
result = "+".join(numbers)

# Print the final, rearranged sum.
print(result)