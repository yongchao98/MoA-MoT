import sys

# Read the entire line from standard input
s = sys.stdin.readline().strip()

# The summands are separated by '+'. We can split the string by '+' to get a list of numbers.
# Since the numbers are single digits '1', '2', '3', '4', sorting them as strings
# is equivalent to sorting them as integers.
summands = s.split('+')

# Sort the list of summands in non-decreasing order.
summands.sort()

# Join the sorted summands back together with '+' and print the result.
# The join method correctly handles printing each number in the final equation.
print('+'.join(summands))