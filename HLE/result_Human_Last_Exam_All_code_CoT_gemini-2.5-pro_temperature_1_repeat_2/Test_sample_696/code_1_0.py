import re

# The user wants to find the 32-character string that replaces "@" in the one-liner.
# The one-liner is: from re import*;print(@r'\d+', input())})or ["NO"])
# The goal is to extract, sort, and print unique numbers from the input string.

# Step 1: The core logic involves finding all numbers. `re.findall(r'\d+', S)` is the way to do this.
# Step 2: The numbers are strings and need to be converted to integers to sort them numerically
# and handle duplicates (e.g., '123' and '0123'). Using `int()`.
# Step 3: We need unique numbers. A `set` is the ideal data structure for this.
# Step 4: The unique numbers must be sorted. The `sorted()` function does this.
# Step 5: If no numbers are found, `re.findall` returns an empty list `[]`. This is a "falsy" value in Python.
# The `... or ["NO"]` structure handles this case. If the left side is falsy, the right side is used.
# Step 6: The output needs to be space-separated numbers, not a Python list.
# `print(*my_list)` unpacks the list and prints its elements with spaces in between.
# `print(*["NO"])` will correctly print `NO`.

# Combining these steps, a functional line of code would be:
# print(*(sorted(set(map(int, re.findall(r'\d+', input())))) or ["NO"]))
# An alternative using a set comprehension, which is often more readable and efficient:
# print(*(sorted({int(n) for n in re.findall(r'\d+', input())}) or ["NO"]))

# Now, we need to fit this logic into the provided template and find the 32 characters for '@'.
# Template: print(@r'\d+', input())})or ["NO"])
# Let's assume the expression is `*sorted({int(n) for n in findall(...)})`.
# The part that should replace '@' is everything before `r'\d+'`.
# This is: `*sorted({int(n) for n in findall(`

# Let's count the characters in this string:
# `*` (1) + `sorted` (6) + `(` (1) + `{` (1) + `int(n)` (5) + ` for ` (5) + `n` (1) + ` in ` (4) + `findall` (7) + `(` (1)
# = 32 characters.
# This matches the requirement.

# The @ replacement string is therefore:
replacement_string = "*sorted({int(n) for n in findall("

print(replacement_string)