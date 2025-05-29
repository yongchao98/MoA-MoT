from itertools import permutations, product

def evaluate_and_format(nums, ops):
    # Helper function to evaluate and format the expression
    def eval_recursive(n, o):
        if len(n) == 1:
            return n[0], str(n[0])
        
        for i in range(len(o)):
            left_val, left_expr = eval_recursive(n[:i+1], o[:i])
            right_val, right_expr = eval_recursive(n[i+1:], o[i+1:])
            
            if o[i] == '+':
                result = left_val + right_val
                expr = f"({left_expr}+{right_expr})"
            elif o[i] == '-':
                result = left_val - right_val
                expr = f"({left_expr}-{right_expr})"
            elif o[i] == '*':
                result = left_val * right_val
                expr = f"({left_expr}*{right_expr})"
            elif o[i] == '/':
                if right_val == 0:
                    continue
                result = left_val / right_val
                expr = f"({left_expr}/{right_expr})"
            
            if result == 24:
                return result, expr
        
        return None, None

    return eval_recursive(nums, ops)

def find_expression_to_24(numbers):
    # Generate all permutations of the numbers
    for num_perm in permutations(numbers):
        # Generate all combinations of operations
        for ops in product('+-*/', repeat=3):
            # Evaluate and format the expression
            result, expression = evaluate_and_format(num_perm, ops)
            if result == 24:
                return f"<<<{expression}>>>"
    return "No solution found"

# Define the numbers
numbers = [1, 1, 2, 11]

# Find and print the expression that results in 24
print(find_expression_to_24(numbers))