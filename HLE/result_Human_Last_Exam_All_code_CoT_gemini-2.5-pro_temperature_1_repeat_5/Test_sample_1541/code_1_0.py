# Read the input sum from the user
s = input()

# Split the string by the '+' delimiter to get a list of number strings
# For example, "3+2+1" becomes ['3', '2', '1']
numbers = s.split('+')

# Sort the list of numbers. Python's default sort works correctly for single-digit strings.
# ['3', '2', '1'] becomes ['1', '2', '3']
numbers.sort()

# Join the sorted list of numbers back into a single string with '+' as the delimiter
# ['1', '2', '3'] becomes "1+2+3"
sorted_sum = "+".join(numbers)

# Print the final, sorted sum. The output is the full equation.
print(sorted_sum)