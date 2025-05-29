from itertools import permutations, product

def evaluate_expression(expr):
    try:
        return eval(expr)
    except ZeroDivisionError:
        return None

def find_expression_to_24(numbers):
    operators = ['+', '-', '*', '/']
    
    # Generate all permutations of the numbers
    for nums in permutations(numbers):
        # Generate all combinations of operators
        for ops in product(operators, repeat=3):
            # Try different groupings with parentheses
            expressions = [
                f"(({nums[0]}{ops[0]}{nums[1]}){ops[1]}{nums[2]}){ops[2]}{nums[3]}",
                f"({nums[0]}{ops[0]}({nums[1]}{ops[1]}{nums[2]})){ops[2]}{nums[3]}",
                f"{nums[0]}{ops[0]}(({nums[1]}{ops[1]}{nums[2]}){ops[2]}{nums[3]})",
                f"{nums[0]}{ops[0]}({nums[1]}{ops[1]}({nums[2]}{ops[2]}{nums[3]}))",
                f"({nums[0]}{ops[0]}{nums[1]}){ops[1]}({nums[2]}{ops[2]}{nums[3]})"
            ]
            
            for expr in expressions:
                if evaluate_expression(expr) == 24:
                    return f"<<<{expr} = 24>>>"
    
    return "No solution found."

# Example usage
numbers = [1, 1, 5, 7]
print(find_expression_to_24(numbers))