from itertools import permutations, product
from sympy import symbols, simplify

def find_expression_to_24(numbers):
    operations = ['+', '-', '*', '/']
    a, b, c, d = symbols('a b c d')
    
    for num_perm in permutations(numbers):
        for ops in product(operations, repeat=3):
            # Generate expressions with different parenthetical groupings
            expressions = [
                ((a, ops[0], b), ops[1], c, ops[2], d),
                (a, ops[0], (b, ops[1], c), ops[2], d),
                (a, ops[0], b, ops[1], (c, ops[2], d)),
                ((a, ops[0], b), ops[1], (c, ops[2], d)),
                (a, ops[0], (b, ops[1], (c, ops[2], d)))
            ]
            for expr in expressions:
                # Construct the expression using sympy
                expr_sympy = simplify(expr)
                # Substitute the numbers into the expression
                expr_sub = expr_sympy.subs({a: num_perm[0], b: num_perm[1], c: num_perm[2], d: num_perm[3]})
                # Evaluate the expression
                if expr_sub == 24:
                    # Format the expression as a string
                    expr_str = str(expr_sympy).replace('a', str(num_perm[0])).replace('b', str(num_perm[1])).replace('c', str(num_perm[2])).replace('d', str(num_perm[3]))
                    return f"<<<{expr_str}>>>"
    return "No solution found"

numbers = [10, 13, 13, 13]
print(find_expression_to_24(numbers))