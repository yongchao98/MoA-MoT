import itertools
import operator

# Define the operations
ops = [operator.add, operator.sub, operator.mul, operator.truediv]
ops_symbols = ['+', '-', '*', '/']

# Function to evaluate an expression with given numbers and operations
def evaluate_expression(nums, ops):
    try:
        # Try different parenthesis placements
        # ((a op b) op c) op d
        result1 = ops[2](ops[1](ops[0](nums[0], nums[1]), nums[2]), nums[3])
        if result1 == 24:
            return f"(({nums[0]}{ops_symbols[0]}{nums[1]}){ops_symbols[1]}{nums[2]}){ops_symbols[2]}{nums[3]}"
        
        # (a op (b op c)) op d
        result2 = ops[2](ops[0](nums[0], ops[1](nums[1], nums[2])), nums[3])
        if result2 == 24:
            return f"({nums[0]}{ops_symbols[0]}({nums[1]}{ops_symbols[1]}{nums[2]})){ops_symbols[2]}{nums[3]}"
        
        # (a op b) op (c op d)
        result3 = ops[1](ops[0](nums[0], nums[1]), ops[2](nums[2], nums[3]))
        if result3 == 24:
            return f"({nums[0]}{ops_symbols[0]}{nums[1]}){ops_symbols[1]}({nums[2]}{ops_symbols[2]}{nums[3]})"
        
        # a op ((b op c) op d)
        result4 = ops[0](nums[0], ops[2](ops[1](nums[1], nums[2]), nums[3]))
        if result4 == 24:
            return f"{nums[0]}{ops_symbols[0]}(({nums[1]}{ops_symbols[1]}{nums[2]}){ops_symbols[2]}{nums[3]})"
        
        # a op (b op (c op d))
        result5 = ops[0](nums[0], ops[1](nums[1], ops[2](nums[2], nums[3])))
        if result5 == 24:
            return f"{nums[0]}{ops_symbols[0]}({nums[1]}{ops_symbols[1]}({nums[2]}{ops_symbols[2]}{nums[3]}))"
        
    except ZeroDivisionError:
        return None

# Input numbers
numbers = [1, 4, 6, 12]

# Try all permutations of numbers and operations
for num_perm in itertools.permutations(numbers):
    for ops_perm in itertools.product(ops, repeat=3):
        result = evaluate_expression(num_perm, ops_perm)
        if result:
            print(f"<<<{result}>>>")
            break