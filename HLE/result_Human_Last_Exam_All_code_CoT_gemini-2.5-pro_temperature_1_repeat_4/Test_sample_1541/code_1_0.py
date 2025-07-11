import sys

# Read the entire line from standard input
s = sys.stdin.readline().strip()

# Split the string by the '+' sign to get a list of numbers
numbers = s.split('+')

# Sort the list of numbers. For single digits, this works correctly.
numbers.sort()

# Join the sorted list back into a string with '+'
# This outputs the final equation with all numbers sorted.
new_sum = '+'.join(numbers)

# Print the final result
print(new_sum)