from itertools import permutations, product

def evaluate_expression(nums, ops):
    # Helper function to evaluate an expression using a stack
    def eval_stack(n, o):
        stack = [n[0]]
        for i in range(3):
            if o[i] == '+':
                stack.append(n[i + 1])
            elif o[i] == '-':
                stack.append(-n[i + 1])
            elif o[i] == '*':
                stack[-1] *= n[i + 1]
            elif o[i] == '/':
                if n[i + 1] == 0:
                    return None
                stack[-1] /= n[i + 1]
        return sum(stack)

    return eval_stack(nums, ops)

def find_expression_to_24(numbers):
    # Generate all permutations of the numbers
    for num_perm in permutations(numbers):
        # Generate all combinations of operations
        for ops in product('+-*/', repeat=3):
            # Evaluate the expression
            if evaluate_expression(num_perm, ops) == 24:
                # Format the expression as a string
                expression = f"({num_perm[0]}{ops[0]}{num_perm[1]}){ops[1]}({num_perm[2]}{ops[2]}{num_perm[3]})"
                return f"<<<{expression}>>>"
    return "No solution found"

# Define the numbers
numbers = [1, 1, 2, 11]

# Find and print the expression that results in 24
print(find_expression_to_24(numbers))