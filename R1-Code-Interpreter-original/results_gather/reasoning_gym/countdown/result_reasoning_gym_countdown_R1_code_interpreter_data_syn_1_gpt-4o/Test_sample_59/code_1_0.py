from itertools import permutations

numbers = [42, 41, 65, 36]
target = 869

# Try all permutations of the numbers
for perm in permutations(numbers):
    a, b, c, d = perm
    # Check different combinations of operations
    if a * b + c + d == target:
        print(f"{a} * {b} + {c} + {d}")
    elif a * b + c - d == target:
        print(f"{a} * {b} + {c} - {d}")
    elif a * b - c + d == target:
        print(f"{a} * {b} - {c} + {d}")
    elif a * b - c - d == target:
        print(f"{a} * {b} - {c} - {d}")
    elif a * b * c + d == target:
        print(f"{a} * {b} * {c} + {d}")
    elif a * b * c - d == target:
        print(f"{a} * {b} * {c} - {d}")
    elif a * b + c * d == target:
        print(f"{a} * {b} + {c} * {d}")
    elif a * b - c * d == target:
        print(f"{a} * {b} - {c} * {d}")