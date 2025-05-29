from itertools import product

# The numbers in the given order
numbers = [2, 2, 7, 2, 3]

# The target result
target = 78

# Possible operations
operations = ['+', '-', '*', '/']

# Function to evaluate an expression
def evaluate_expression(numbers, ops):
    expression = f"{numbers[0]}"
    for i in range(1, len(numbers)):
        expression += f"{ops[i-1]}{numbers[i]}"
    try:
        return eval(expression)
    except ZeroDivisionError:
        return None

# Try all combinations of operations
for ops in product(operations, repeat=4):
    if evaluate_expression(numbers, ops) == target:
        result = [numbers[0]]
        for i in range(4):
            result.append(ops[i])
            result.append(numbers[i+1])
        print(result)
        break