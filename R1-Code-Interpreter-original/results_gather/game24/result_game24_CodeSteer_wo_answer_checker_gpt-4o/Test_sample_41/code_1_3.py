from itertools import permutations, product
from sympy import symbols, Eq, solve

def find_expression_to_24(numbers):
    operations = ['+', '-', '*', '/']
    for num_perm in permutations(numbers):
        for ops in product(operations, repeat=3):
            # Define symbols for the numbers
            a, b, c, d = symbols('a b c d')
            # Generate expressions with different parenthetical groupings
            expressions = [
                ((a, ops[0], b), ops[1], c, ops[2], d),
                (a, ops[0], (b, ops[1], c), ops[2], d),
                (a, ops[0], b, ops[1], (c, ops[2], d)),
                ((a, ops[0], b), ops[1], (c, ops[2], d)),
                (a, ops[0], (b, ops[1], (c, ops[2], d)))
            ]
            for expr in expressions:
                # Convert the tuple expression to a string
                expr_str = str(expr).replace(',', '').replace("'", '')
                # Substitute the numbers into the expression
                expr_eval = expr_str.replace('a', str(num_perm[0])).replace('b', str(num_perm[1])).replace('c', str(num_perm[2])).replace('d', str(num_perm[3]))
                # Evaluate the expression
                try:
                    if eval(expr_eval) == 24:
                        return f"<<<{expr_eval}>>>"
                except ZeroDivisionError:
                    continue
    return "No solution found"

numbers = [10, 13, 13, 13]
print(find_expression_to_24(numbers))