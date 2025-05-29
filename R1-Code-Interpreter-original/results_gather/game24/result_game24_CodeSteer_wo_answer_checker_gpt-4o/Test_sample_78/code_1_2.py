import itertools
import operator

# Define the numbers and target
numbers = [3, 6, 6, 12]
target = 24

# Define possible operations and their symbols
operations = [operator.add, operator.sub, operator.mul, operator.truediv]
operation_symbols = {operator.add: '+', operator.sub: '-', operator.mul: '*', operator.truediv: '/'}

def evaluate_expression(nums, ops):
    # Try different parenthesizations
    try:
        # ((a op1 b) op2 c) op3 d
        if ops[2](ops[1](ops[0](nums[0], nums[1]), nums[2]), nums[3]) == target:
            return f"<<<(({nums[0]}{operation_symbols[ops[0]]}{nums[1]}){operation_symbols[ops[1]]}{nums[2]}){operation_symbols[ops[2]]}{nums[3]}={target}>>>"
        # (a op1 (b op2 c)) op3 d
        if ops[2](ops[0](nums[0], ops[1](nums[1], nums[2])), nums[3]) == target:
            return f"<<<({nums[0]}{operation_symbols[ops[0]]}({nums[1]}{operation_symbols[ops[1]]}{nums[2]})){operation_symbols[ops[2]]}{nums[3]}={target}>>>"
        # (a op1 b) op2 (c op3 d)
        if ops[1](ops[0](nums[0], nums[1]), ops[2](nums[2], nums[3])) == target:
            return f"<<<({nums[0]}{operation_symbols[ops[0]]}{nums[1]}){operation_symbols[ops[1]]}({nums[2]}{operation_symbols[ops[2]]}{nums[3]})={target}>>>"
        # a op1 ((b op2 c) op3 d)
        if ops[0](nums[0], ops[2](ops[1](nums[1], nums[2]), nums[3])) == target:
            return f"<<<{nums[0]}{operation_symbols[ops[0]]}(({nums[1]}{operation_symbols[ops[1]]}{nums[2]}){operation_symbols[ops[2]]}{nums[3]})={target}>>>"
        # a op1 (b op2 (c op3 d))
        if ops[0](nums[0], ops[1](nums[1], ops[2](nums[2], nums[3]))) == target:
            return f"<<<{nums[0]}{operation_symbols[ops[0]]}({nums[1]}{operation_symbols[ops[1]]}({nums[2]}{operation_symbols[ops[2]]}{nums[3]}))={target}>>>"
    except ZeroDivisionError:
        pass
    return None

# Try all permutations of numbers and operations
for num_perm in itertools.permutations(numbers):
    for ops in itertools.product(operations, repeat=3):
        result = evaluate_expression(num_perm, ops)
        if result:
            print(result)
            break