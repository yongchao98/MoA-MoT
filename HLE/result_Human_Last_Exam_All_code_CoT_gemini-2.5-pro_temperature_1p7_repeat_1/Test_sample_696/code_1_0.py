from re import *

# The user wants to find the 32 characters to replace '@' in the code:
# from re import*;print(@r'\d+', input())})or ["NO"])

# Let's break down the logic needed:
# 1. Find all substrings that are numbers. `findall(r'\d+', input())` does this.
# 2. Convert these number strings to integers to handle things like leading zeros ('07' -> 7) and allow numerical sorting.
# 3. Get only the unique numbers. A `set` is perfect for this.
# 4. Sort the unique numbers. The `sorted()` function does this.
# 5. A set comprehension `{int(i) for i in ...}` is a concise way to handle steps 2 and 3.
# 6. So, the core logic is `sorted({int(i) for i in findall(r'\d+', input())})`. This produces a sorted list of unique integers.
# 7. To handle the "NO" case, we can use the `or` operator. If the list of numbers is empty (which is falsy), `... or ["NO"]` will yield `["NO"]`.
# 8. To print the numbers separated by spaces, we use the `*` unpacking operator: `print(*result)`.

# The original template `print(@r'\d+', input())})or ["NO"])` has a syntax error due to a mismatched `}`.
# Assuming this is a typo and constructing a valid one-liner with the above logic, we get:
# print(*(sorted({int(i)for i in findall(r'\d+', input())}) or ["NO"]))

# The part of this code that starts the logic and uses `findall` corresponds to `@` in the prompt.
# The string `sorted({int(i)for i in findall` is exactly 32 characters.
# This appears to be the intended answer to the puzzle.

# The following code is the complete, corrected one-liner that solves the problem.
# The output of this script will be the sorted unique numbers or "NO".

print(*(sorted({int(i)for i in findall(r'\d+', input())}) or ["NO"]))