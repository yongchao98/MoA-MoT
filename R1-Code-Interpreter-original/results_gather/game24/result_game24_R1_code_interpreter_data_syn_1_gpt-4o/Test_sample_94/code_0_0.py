from itertools import permutations, product

def evaluate_expression(nums, ops):
    # Try different parenthesis placements
    try:
        # ((a op1 b) op2 c) op3 d
        result1 = ops[2](ops[1](ops[0](nums[0], nums[1]), nums[2]), nums[3])
        if result1 == 24:
            return f"(({nums[0]}{ops[0].__name__}{nums[1]}){ops[1].__name__}{nums[2]}){ops[2].__name__}{nums[3]}"
        
        # (a op1 (b op2 c)) op3 d
        result2 = ops[2](ops[0](nums[0], ops[1](nums[1], nums[2])), nums[3])
        if result2 == 24:
            return f"({nums[0]}{ops[0].__name__}({nums[1]}{ops[1].__name__}{nums[2]})){ops[2].__name__}{nums[3]}"
        
        # (a op1 b) op2 (c op3 d)
        result3 = ops[1](ops[0](nums[0], nums[1]), ops[2](nums[2], nums[3]))
        if result3 == 24:
            return f"({nums[0]}{ops[0].__name__}{nums[1]}){ops[1].__name__}({nums[2]}{ops[2].__name__}{nums[3]})"
        
        # a op1 ((b op2 c) op3 d)
        result4 = ops[0](nums[0], ops[2](ops[1](nums[1], nums[2]), nums[3]))
        if result4 == 24:
            return f"{nums[0]}{ops[0].__name__}(({nums[1]}{ops[1].__name__}{nums[2]}){ops[2].__name__}{nums[3]})"
        
        # a op1 (b op2 (c op3 d))
        result5 = ops[0](nums[0], ops[1](nums[1], ops[2](nums[2], nums[3])))
        if result5 == 24:
            return f"{nums[0]}{ops[0].__name__}({nums[1]}{ops[1].__name__}({nums[2]}{ops[2].__name__}{nums[3]}))"
        
    except ZeroDivisionError:
        return None

    return None

# Define operations
operations = [lambda x, y: x + y, lambda x, y: x - y, lambda x, y: x * y, lambda x, y: x / y]

# Map operation names for output
for op in operations:
    op.__name__ = {lambda x, y: x + y: '+', lambda x, y: x - y: '-', lambda x, y: x * y: '*', lambda x, y: x / y: '/'}[op]

# Numbers to use
numbers = [1, 3, 7, 8]

# Try all permutations of numbers and operations
for num_perm in permutations(numbers):
    for ops in product(operations, repeat=3):
        expression = evaluate_expression(num_perm, ops)
        if expression:
            print(expression)
            break