# Plan: The one-liner needs to find all numbers, make them unique, convert to integers, sort them,
# and then print them space-separated. If no numbers are found, it should print "NO".
#
# 1. `findall(r'\d+', input())` finds all number strings.
# 2. A set comprehension `{int(n) for n in ...}` converts them to unique integers.
# 3. `sorted()` sorts the resulting set.
# 4. The `*` operator unpacks the list for `print`.
# 5. The `or ["NO"]` handles the case with no numbers.
#
# The code that must precede `r'\d+', input())})or ["NO"])` to make this work is
# `*sorted({int(n) for n in findall(`.
# Let's count the characters:
# '*': 1
# 'sorted': 6
# '(': 1
# '{': 1
# 'int(n)': 5
# ' ': 1
# 'for': 3
# ' ': 1
# 'n': 1
# ' ': 1
# 'in': 2
# ' ': 1
# 'findall': 7
# '(': 1
# Total: 1+6+1+1+5+1+3+1+1+1+2+1+7+1 = 32 characters.

solution_string = "*sorted({int(n) for n in findall("
print(solution_string)