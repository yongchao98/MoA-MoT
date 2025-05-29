from sympy import symbols, simplify
from itertools import permutations, product

def find_expression(numbers, target):
    if len(numbers) == 1:
        if simplify(numbers[0]) == target:
            return str(numbers[0])
        else:
            return None

    for i, j in permutations(range(len(numbers)), 2):
        num1, num2 = numbers[i], numbers[j]
        remaining = [numbers[k] for k in range(len(numbers)) if k != i and k != j]

        for op in ['+', '-', '*', '/']:
            if op == '/' and num2 == 0:
                continue

            expr = f"({num1}{op}{num2})"
            new_numbers = [expr] + remaining
            result = find_expression(new_numbers, target)

            if result:
                return result

    return None

numbers = [3, 5, 5, 12]
target = 24
x = symbols('x')
expression = find_expression([str(x) for x in numbers], target)

if expression:
    print(f"<<<{expression}>>>")
else:
    print("No solution found.")