from itertools import product

def evaluate_expression_with_parentheses(numbers, ops, parens):
    expression = ""
    num_index = 0
    op_index = 0
    for i in range(len(parens)):
        if parens[i] == '(':
            expression += '('
        expression += str(numbers[num_index])
        num_index += 1
        if op_index < len(ops):
            expression += ops[op_index]
            op_index += 1
        if parens[i] == ')':
            expression += ')'
    expression += str(numbers[num_index])  # Add the last number
    try:
        return eval(expression)
    except ZeroDivisionError:
        return None

def find_expression_with_parentheses(numbers, target):
    operators = ['+', '-', '*', '/']
    n = len(numbers)
    for ops in product(operators, repeat=n - 1):
        for parens in product(['', '(', ')'], repeat=n - 1):
            if parens.count('(') == parens.count(')'):
                result = evaluate_expression_with_parentheses(numbers, ops, parens)
                if result == target:
                    return [numbers[0]] + [item for pair in zip(ops, numbers[1:]) for item in pair]
    return None

numbers = [6, 9, 6, 7, 2, 8]
target = 100
result = find_expression_with_parentheses(numbers, target)
print(result)