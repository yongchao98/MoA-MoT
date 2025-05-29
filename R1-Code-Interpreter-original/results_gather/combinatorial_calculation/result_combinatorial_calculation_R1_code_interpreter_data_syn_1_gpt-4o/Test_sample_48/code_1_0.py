from itertools import product

def evaluate_expression(numbers, ops, parens):
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
    except (ZeroDivisionError, SyntaxError):
        return None

def find_expression(numbers, target):
    operators = ['+', '-', '*', '/']
    n = len(numbers)
    for ops in product(operators, repeat=n - 1):
        for parens in product(['', '(', ')'], repeat=n - 1):
            if parens.count('(') == parens.count(')'):
                # Ensure no operator is directly before a closing parenthesis or after an opening parenthesis
                valid = True
                for j in range(len(parens) - 1):
                    if parens[j] == '(' and ops[j] in operators:
                        valid = False
                        break
                    if parens[j] == ')' and ops[j] in operators:
                        valid = False
                        break
                if not valid:
                    continue
                result = evaluate_expression(numbers, ops, parens)
                if result == target:
                    return [numbers[0]] + [item for pair in zip(ops, numbers[1:]) for item in pair]
    return None

numbers = [6, 9, 6, 7, 2, 8]
target = 100
result = find_expression(numbers, target)
print(result)