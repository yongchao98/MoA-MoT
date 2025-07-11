from re import findall

# Read the input string from the user.
input_string = input()

# Find all sequences of digits.
numbers_as_strings = findall(r'\d+', input_string)

# If numbers were found, process them. Otherwise, prepare the "NO" message.
# The 'or ["NO"]' trick handles the case where the initial list is empty.
# 1. Create a set to get unique number strings: set(numbers_as_strings)
# 2. Sort the set numerically: sorted(..., key=int)
# 3. The result is either a list of sorted number strings or ["NO"].
# 4. Join the elements of the list with spaces.
output_list = sorted(set(numbers_as_strings), key=int) or ["NO"]
print(' '.join(output_list))
