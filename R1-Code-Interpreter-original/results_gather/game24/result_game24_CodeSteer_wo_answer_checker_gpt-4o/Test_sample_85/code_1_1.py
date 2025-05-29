import itertools

def find_expression_to_24(numbers):
    operators = ['+', '-', '*', '/']
    
    # Generate all permutations of the numbers
    for num_perm in itertools.permutations(numbers):
        # Generate all combinations of operators
        for ops in itertools.product(operators, repeat=3):
            # Generate all possible parenthesized expressions
            expressions = [
                f"(({num_perm[0]}{ops[0]}{num_perm[1]}){ops[1]}{num_perm[2]}){ops[2]}{num_perm[3]}",
                f"({num_perm[0]}{ops[0]}({num_perm[1]}{ops[1]}{num_perm[2]})){ops[2]}{num_perm[3]}",
                f"{num_perm[0]}{ops[0]}(({num_perm[1]}{ops[1]}{num_perm[2]}){ops[2]}{num_perm[3]})",
                f"{num_perm[0]}{ops[0]}({num_perm[1]}{ops[1]}({num_perm[2]}{ops[2]}{num_perm[3]}))",
                f"({num_perm[0]}{ops[0]}{num_perm[1]}){ops[1]}({num_perm[2]}{ops[2]}{num_perm[3]})"
            ]
            
            for expr in expressions:
                try:
                    # Evaluate the expression and check if it equals 24
                    if abs(eval(expr) - 24) < 1e-9:  # Use a tolerance for floating-point arithmetic
                        return f"<<<{expr}>>>"
                except ZeroDivisionError:
                    continue

# Input numbers
numbers = [1, 7, 10, 12]
result = find_expression_to_24(numbers)
print(result)