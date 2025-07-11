def solve_equation():
    """
    This function solves the final decoded equation.
    The puzzle's final word 'anhaeketn' decodes to the numbers 8 and 3.
    Assuming the operation is multiplication.
    """
    number1 = 8
    number2 = 3
    operator = "*"
    result = number1 * number2
    
    print(f"{number1} {operator} {number2} = {result}")

solve_equation()