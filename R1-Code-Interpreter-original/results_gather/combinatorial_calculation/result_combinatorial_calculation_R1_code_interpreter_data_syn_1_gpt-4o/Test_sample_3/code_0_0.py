from itertools import product

numbers = [5, 2, 4]
target = 10
operations = ['+', '-', '*', '/']

def evaluate_expression(expr):
    try:
        return eval(expr)
    except ZeroDivisionError:
        return None

def find_expression(numbers, target):
    for ops in product(operations, repeat=2):
        # Without parentheses
        expr1 = f"{numbers[0]} {ops[0]} {numbers[1]} {ops[1]} {numbers[2]}"
        if evaluate_expression(expr1) == target:
            return expr1

        # With parentheses
        expr2 = f"({numbers[0]} {ops[0]} {numbers[1]}) {ops[1]} {numbers[2]}"
        if evaluate_expression(expr2) == target:
            return expr2

        expr3 = f"{numbers[0]} {ops[0]} ({numbers[1]} {ops[1]} {numbers[2]})"
        if evaluate_expression(expr3) == target:
            return expr3

    return None

expression = find_expression(numbers, target)
print(expression)