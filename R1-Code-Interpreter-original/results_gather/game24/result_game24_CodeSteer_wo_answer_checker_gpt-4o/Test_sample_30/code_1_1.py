import itertools
import operator

# Define the numbers and target
numbers = [6, 6, 10, 13]
target = 24

# Define possible operations
ops = [operator.add, operator.sub, operator.mul, operator.truediv]
op_symbols = ['+', '-', '*', '/']

def evaluate_expression(nums, ops):
    """Evaluate the expression with given numbers and operations."""
    try:
        # Try different groupings with parentheses
        result1 = ops[0](nums[0], ops[1](nums[1], ops[2](nums[2], nums[3])))
        result2 = ops[0](ops[1](nums[0], nums[1]), ops[2](nums[2], nums[3]))
        result3 = ops[0](ops[1](ops[2](nums[0], nums[1]), nums[2]), nums[3])
        result4 = ops[0](nums[0], ops[1](ops[2](nums[1], nums[2]), nums[3]))
        result5 = ops[0](ops[1](nums[0], ops[2](nums[1], nums[2])), nums[3])
        
        return [result1, result2, result3, result4, result5]
    except ZeroDivisionError:
        return [None] * 5

def find_expression():
    """Find an expression that evaluates to the target."""
    for num_perm in itertools.permutations(numbers):
        for op_perm in itertools.product(ops, repeat=3):
            results = evaluate_expression(num_perm, op_perm)
            for i, result in enumerate(results):
                if result == target:
                    # Format the expression based on the grouping
                    if i == 0:
                        expr = f"({num_perm[0]} {op_symbols[ops.index(op_perm[0])]} ({num_perm[1]} {op_symbols[ops.index(op_perm[1])]} ({num_perm[2]} {op_symbols[ops.index(op_perm[2])]} {num_perm[3]})))"
                    elif i == 1:
                        expr = f"(({num_perm[0]} {op_symbols[ops.index(op_perm[0])]} {num_perm[1]}) {op_symbols[ops.index(op_perm[1])]} ({num_perm[2]} {op_symbols[ops.index(op_perm[2])]} {num_perm[3]}))"
                    elif i == 2:
                        expr = f"(({num_perm[0]} {op_symbols[ops.index(op_perm[0])]} ({num_perm[1]} {op_symbols[ops.index(op_perm[1])]} {num_perm[2]})) {op_symbols[ops.index(op_perm[2])]} {num_perm[3]})"
                    elif i == 3:
                        expr = f"({num_perm[0]} {op_symbols[ops.index(op_perm[0])]} (({num_perm[1]} {op_symbols[ops.index(op_perm[1])]} {num_perm[2]}) {op_symbols[ops.index(op_perm[2])]} {num_perm[3]}))"
                    elif i == 4:
                        expr = f"(({num_perm[0]} {op_symbols[ops.index(op_perm[0])]} {num_perm[1]}) {op_symbols[ops.index(op_perm[1])]} ({num_perm[2]} {op_symbols[ops.index(op_perm[2])]} {num_perm[3]}))"
                    
                    return f"<<<{expr} = {target}>>>"
    return "No solution found."

# Find and print the expression
print(find_expression())