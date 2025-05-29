from itertools import permutations

numbers = [82, 38, 56, 49]
target = 588

def find_expression(numbers, target):
    for perm in permutations(numbers):
        a, b, c, d = perm
        # Check different combinations of operations
        if a * b + c + d == target:
            return f"{a} * {b} + {c} + {d}"
        if a * b + c - d == target:
            return f"{a} * {b} + {c} - {d}"
        if a * b - c + d == target:
            return f"{a} * {b} - {c} + {d}"
        if a * b - c - d == target:
            return f"{a} * {b} - {c} - {d}"
        if a * b * c + d == target:
            return f"{a} * {b} * {c} + {d}"
        if a * b * c - d == target:
            return f"{a} * {b} * {c} - {d}"
        if a * b * c * d == target:
            return f"{a} * {b} * {c} * {d}"
        if a * b * c / d == target:
            return f"{a} * {b} * {c} / {d}"
        if a * b / c + d == target:
            return f"{a} * {b} / {c} + {d}"
        if a * b / c - d == target:
            return f"{a} * {b} / {c} - {d}"
        if a / b * c + d == target:
            return f"{a} / {b} * {c} + {d}"
        if a / b * c - d == target:
            return f"{a} / {b} * {c} - {d}"
        if a / b / c * d == target:
            return f"{a} / {b} / {c} * {d}"
        if a / b / c / d == target:
            return f"{a} / {b} / {c} / {d}"
    return None

expression = find_expression(numbers, target)
print(expression)