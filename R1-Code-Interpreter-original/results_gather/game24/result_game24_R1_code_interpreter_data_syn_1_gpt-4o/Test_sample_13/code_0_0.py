from itertools import permutations
import operator

# Define the operations
ops = {
    '+': operator.add,
    '-': operator.sub,
    '*': operator.mul,
    '/': operator.truediv
}

# Function to evaluate an expression
def eval_expr(expr):
    try:
        return eval(expr)
    except ZeroDivisionError:
        return None

# Generate all permutations of the numbers
numbers = [3, 8, 10, 10]
perms = permutations(numbers)

# Generate all possible expressions
for perm in perms:
    a, b, c, d = perm
    for op1 in ops:
        for op2 in ops:
            for op3 in ops:
                # Form different expressions
                expr1 = f"({a}{op1}{b}){op2}({c}{op3}{d})"
                expr2 = f"(({a}{op1}{b}){op2}{c}){op3}{d}"
                expr3 = f"({a}{op1}({b}{op2}{c})){op3}{d}"
                expr4 = f"{a}{op1}(({b}{op2}{c}){op3}{d})"
                expr5 = f"{a}{op1}({b}{op2}({c}{op3}{d}))"
                
                # Evaluate expressions
                for expr in [expr1, expr2, expr3, expr4, expr5]:
                    if eval_expr(expr) == 24:
                        print(expr)
                        break