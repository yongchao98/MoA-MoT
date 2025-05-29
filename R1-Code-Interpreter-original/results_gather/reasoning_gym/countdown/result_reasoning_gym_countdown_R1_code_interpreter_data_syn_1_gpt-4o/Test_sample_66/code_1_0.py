# Let's try different combinations to make 16 using 54, 46, and 6
from itertools import permutations

numbers = [54, 46, 6]
target = 16

# Try all permutations of the numbers with different operations
def evaluate_expression(a, b, c):
    # Try different operations to see if we can reach the target
    if a + b - c == target:
        return f"{a} + {b} - {c}"
    if a - b + c == target:
        return f"{a} - {b} + {c}"
    if a * b - c == target:
        return f"{a} * {b} - {c}"
    if a - b * c == target:
        return f"{a} - {b} * {c}"
    if a + b + c == target:
        return f"{a} + {b} + {c}"
    if a * b + c == target:
        return f"{a} * {b} + {c}"
    if a + b * c == target:
        return f"{a} + {b} * {c}"
    if a * b * c == target:
        return f"{a} * {b} * {c}"
    return None

# Check all permutations
for perm in permutations(numbers):
    result = evaluate_expression(*perm)
    if result:
        print(result)
        break