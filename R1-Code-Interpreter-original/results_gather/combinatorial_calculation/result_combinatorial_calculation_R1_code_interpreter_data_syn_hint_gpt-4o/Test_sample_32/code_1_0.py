from sympy import symbols, Add, Mul, Rational

# Define the numbers as symbols
a, b, c, d, e = symbols('a b c d e')

# Given numbers
numbers = [1, 9, 4, 9, 7]

# Possible operations
operations = [Add, Mul, Rational]

# Function to evaluate expression
def evaluate_expression(expr):
    try:
        return expr.evalf()
    except ZeroDivisionError:
        return None

# Generate expressions with operations and parentheses
def generate_expressions(numbers, ops):
    a, b, c, d, e = numbers
    for op1 in ops:
        for op2 in ops:
            for op3 in ops:
                for op4 in ops:
                    expr1 = op1(a, op2(b, op3(c, op4(d, e))))
                    expr2 = op1(op2(a, b), op3(c, op4(d, e)))
                    expr3 = op1(op2(op3(a, b), c), op4(d, e))
                    expr4 = op1(op2(a, op3(b, c)), op4(d, e))
                    expr5 = op1(op2(a, b), op3(op4(c, d), e))
                    for expr in [expr1, expr2, expr3, expr4, expr5]:
                        if evaluate_expression(expr) == 69:
                            return expr
    return None

# Try all combinations of operations
result = generate_expressions(numbers, operations)
if result:
    print(result)