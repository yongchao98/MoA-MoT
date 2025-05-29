from itertools import product

def evaluate_expression(numbers, operations):
    expression = str(numbers[0])
    for i in range(len(operations)):
        expression += operations[i] + str(numbers[i + 1])
    try:
        return eval(expression)
    except ZeroDivisionError:
        return None

def find_expression(numbers, target):
    operators = ['+', '-', '*', '/']
    for ops in product(operators, repeat=len(numbers) - 1):
        if evaluate_expression(numbers, ops) == target:
            return [numbers[0]] + [item for pair in zip(ops, numbers[1:]) for item in pair]
    return None

numbers = [6, 9, 6, 7, 2, 8]
target = 100
result = find_expression(numbers, target)
print(result)