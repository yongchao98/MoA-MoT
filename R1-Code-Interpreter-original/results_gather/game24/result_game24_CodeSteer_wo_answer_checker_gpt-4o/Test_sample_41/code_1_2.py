from itertools import permutations, product

def evaluate_expression(expr):
    try:
        return eval(expr)
    except ZeroDivisionError:
        return None

def find_expression_to_24(numbers):
    operations = ['+', '-', '*', '/']
    for num_perm in permutations(numbers):
        for ops in product(operations, repeat=3):
            # Generate expressions with different parenthetical groupings
            expressions = [
                f"(({num_perm[0]}{ops[0]}{num_perm[1]}){ops[1]}{num_perm[2]}){ops[2]}{num_perm[3]}",
                f"({num_perm[0]}{ops[0]}({num_perm[1]}{ops[1]}{num_perm[2]})){ops[2]}{num_perm[3]}",
                f"{num_perm[0]}{ops[0]}(({num_perm[1]}{ops[1]}{num_perm[2]}){ops[2]}{num_perm[3]})",
                f"{num_perm[0]}{ops[0]}({num_perm[1]}{ops[1]}({num_perm[2]}{ops[2]}{num_perm[3]}))",
                f"({num_perm[0]}{ops[0]}{num_perm[1]}){ops[1]}({num_perm[2]}{ops[2]}{num_perm[3]})"
            ]
            for expr in expressions:
                if evaluate_expression(expr) == 24:
                    return f"<<<{expr}>>>"
    return "No solution found"

numbers = [10, 13, 13, 13]
print(find_expression_to_24(numbers))