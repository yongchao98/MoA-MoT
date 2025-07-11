import sys

# Read the sum from the input line.
# Example input: "3+2+1+4"
input_sum = sys.stdin.readline().strip()

# The numbers are separated by '+'. We can split the string by '+'
# to get a list of the number strings.
# Example: ['3', '2', '1', '4']
numbers = input_sum.split('+')

# Sort the list of numbers. Since they are single digits,
# a standard string sort will correctly order them non-decreasingly.
# Example: ['1', '2', '3', '4']
numbers.sort()

# Join the sorted list of numbers back into a single string,
# with '+' as the separator between them.
# Example: "1+2+3+4"
sorted_sum = '+'.join(numbers)

# Print the final rearranged sum. Each number and the '+' signs
# are part of the final equation string.
print(sorted_sum)