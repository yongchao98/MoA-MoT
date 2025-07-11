# Final plan:
# 1. Use findall from the re module to extract all digit sequences.
# 2. Use a set comprehension {int(i) for i in ...} to convert them to integers and get unique values.
# 3. Sort the resulting set to get an ordered list of numbers.
# 4. Use the `*` operator to print the list elements space-separated.
# 5. Use the `or ["NO"]` idiom to provide a default value if the list of numbers is empty.

# The 32-character code snippet will be constructed to initiate this sequence of operations.
# The code to be filled in is: *sorted({int(i) for i in findall

# Let's print the reconstructed full line of code to demonstrate.
# Note: The original prompt's one-liner has syntactic issues. 
# We present the corrected, functional version that our 32-character snippet builds towards.
s = input()
from re import findall

# This is the list of numbers, or an empty list if none are found
numbers = sorted({int(i) for i in findall(r'\d+', s)})

# The 'or' provides a default. The '*' unpacks the list for printing.
# If numbers = [123, 456], it prints 123 456
# If numbers = [], it becomes *["NO"], which prints NO
print(*(numbers or ["NO"]))