import itertools
import operator

# Define the numbers and target
numbers = [8, 9, 10, 13]
target = 24

# Define possible operations
operations = [operator.add, operator.sub, operator.mul, operator.truediv]
operation_symbols = ['+', '-', '*', '/']

def evaluate_expression(nums, ops):
    """Evaluate the expression with given numbers and operations."""
    try:
        # Try different parenthesizations
        result1 = ops[0](ops[1](nums[0], nums[1]), ops[2](nums[2], nums[3]))
        result2 = ops[0](nums[0], ops[1](nums[1], ops[2](nums[2], nums[3])))
        result3 = ops[0](ops[1](ops[2](nums[0], nums[1]), nums[2]), nums[3])
        result4 = ops[0](nums[0], ops[1](ops[2](nums[1], nums[2]), nums[3]))
        result5 = ops[0](ops[1](nums[0], ops[2](nums[1], nums[2])), nums[3])
        
        return [result1, result2, result3, result4, result5]
    except ZeroDivisionError:
        return [None] * 5

def find_solution():
    for num_perm in itertools.permutations(numbers):
        for ops in itertools.product(operations, repeat=3):
            results = evaluate_expression(num_perm, ops)
            for i, result in enumerate(results):
                if result == target:
                    # Format the solution
                    if i == 0:
                        solution = f"(({num_perm[0]} {operation_symbols[ops[0]]} {num_perm[1]}) {operation_symbols[ops[1]]} ({num_perm[2]} {operation_symbols[ops[2]]} {num_perm[3]}))"
                    elif i == 1:
                        solution = f"({num_perm[0]} {operation_symbols[ops[0]]} ({num_perm[1]} {operation_symbols[ops[1]]} ({num_perm[2]} {operation_symbols[ops[2]]} {num_perm[3]})))"
                    elif i == 2:
                        solution = f"(({num_perm[0]} {operation_symbols[ops[0]]} ({num_perm[1]} {operation_symbols[ops[1]]} {num_perm[2]})) {operation_symbols[ops[2]]} {num_perm[3]})"
                    elif i == 3:
                        solution = f"({num_perm[0]} {operation_symbols[ops[0]]} ({num_perm[1]} {operation_symbols[ops[1]]} ({num_perm[2]} {operation_symbols[ops[2]]} {num_perm[3]})))"
                    elif i == 4:
                        solution = f"(({num_perm[0]} {operation_symbols[ops[0]]} {num_perm[1]}) {operation_symbols[ops[1]]} ({num_perm[2]} {operation_symbols[ops[2]]} {num_perm[3]}))"
                    
                    print(f"<<<{solution} = {target}>>>")
                    return

find_solution()