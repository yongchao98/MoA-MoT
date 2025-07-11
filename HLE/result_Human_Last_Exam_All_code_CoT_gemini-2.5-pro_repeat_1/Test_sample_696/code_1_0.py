# Final working code for the problem.
# The 32-character string for @ is ` *( sorted(set(map(int, findall(`
# The full, corrected one-liner is constructed below.

from re import *

# The problem asks for the 32 characters that replace '@'.
# Let's break down the logic:
# 1. findall(r'\d+', input()): Finds all number strings.
#    e.g., for "hello 123 world 456 hello 123", this gives ['123', '456', '123']
# 2. map(int, ...): Converts them to integers.
#    e.g., map object for 123, 456, 123
# 3. set(...): Gets the unique numbers.
#    e.g., {123, 456}
# 4. sorted(...): Sorts them.
#    e.g., [123, 456]
# 5. ... or ["NO"]: If the list is empty, it becomes ["NO"].
# 6. print(*...): Prints the list elements separated by spaces.
#    e.g., print(123, 456) -> "123 456"
#    e.g., print("NO") -> "NO"

# The original one-liner has syntax errors. Here is the corrected, functional version.
# The expression that replaces '@' and the subsequent parentheses/braces is built from the logic above.
print(*( sorted(set(map(int, findall(r'\d+', input())))) or ["NO"]))
