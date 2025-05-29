from itertools import permutations, product

def evaluate_expression(nums, ops):
    # Recursive function to evaluate expressions with different groupings
    def eval_recursive(n, o):
        if len(n) == 1:
            return n[0]
        for i in range(len(o)):
            left = eval_recursive(n[:i+1], o[:i])
            right = eval_recursive(n[i+1:], o[i+1:])
            if o[i] == '+':
                result = left + right
            elif o[i] == '-':
                result = left - right
            elif o[i] == '*':
                result = left * right
            elif o[i] == '/':
                if right == 0:
                    continue
                result = left / right
            if result == 24:
                return result
        return None

    return eval_recursive(nums, ops)

def find_expression_to_24(numbers):
    # Generate all permutations of the numbers
    for num_perm in permutations(numbers):
        # Generate all combinations of operations
        for ops in product('+-*/', repeat=3):
            # Evaluate the expression
            if evaluate_expression(num_perm, ops) == 24:
                # Format the expression as a string
                expression = f"(({num_perm[0]}{ops[0]}{num_perm[1]}){ops[1]}{num_perm[2]}){ops[2]}{num_perm[3]}"
                return f"<<<{expression}>>>"
    return "No solution found"

# Define the numbers
numbers = [1, 1, 2, 11]

# Find and print the expression that results in 24
print(find_expression_to_24(numbers))