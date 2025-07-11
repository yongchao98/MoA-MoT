from re import *;print(*(sorted({int(n)for n in findall(r'\d+', input())}))or ["NO"]))
# The 32 characters that replace the '@' are: *(sorted({int(n)for n in findall(
# Here is the breakdown of the logic for the final output:
# 1. findall(r'\d+', input()) finds all number strings. e.g., for "100 apples 50", it's ['100', '50']
# 2. The set comprehension {int(n) for n in ...} converts them to unique integers. e.g., {100, 50}
# 3. sorted(...) sorts the set into a list. e.g., [50, 100]
# 4. If the list is empty (no numbers found), it's 'falsy', so the `or ["NO"]` part makes the expression result in ["NO"].
# 5. The '*' operator unpacks the list elements as arguments to print.
#    - print(*[50, 100]) becomes print(50, 100), which outputs "50 100"
#    - print(*["NO"]) becomes print("NO"), which outputs "NO"