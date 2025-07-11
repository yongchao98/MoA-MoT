# This script solves the lexicographical string comparison problem.

# The memory calculation for the most efficient C program is as follows:
# To store two strings of maximum length 100, we need two char arrays.
# Each array needs 100 chars + 1 null terminator = 101 bytes.
# We also need a loop counter. Since the length is at most 100, a 'char' type (1 byte) is sufficient.
# Final Equation: m = memory(string1) + memory(string2) + memory(loop_counter)
# The numbers in the equation are:
# m = 101 + 101 + 1
# m = 203 bytes

# Read the two strings from the user
s1 = input()
s2 = input()

# Convert both strings to lowercase for case-insensitive comparison
s1_lower = s1.lower()
s2_lower = s2.lower()

# Perform the lexicographical comparison and print the result
if s1_lower < s2_lower:
    print("-1")
elif s1_lower > s2_lower:
    print("1")
else:
    print("0")
