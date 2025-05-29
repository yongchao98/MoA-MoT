from itertools import permutations
import operator

# Define the operations
ops = [('+', operator.add), ('-', operator.sub), ('*', operator.mul), ('/', operator.truediv)]

def evaluate_expression(nums, ops):
    """Evaluate the expression with the given numbers and operations."""
    try:
        # Try different groupings by evaluating the expression in different orders
        result1 = ops[1][1](ops[0][1](nums[0], nums[1]), ops[2][1](nums[2], nums[3]))
        result2 = ops[0][1](nums[0], ops[1][1](nums[1], ops[2][1](nums[2], nums[3])))
        result3 = ops[0][1](ops[1][1](nums[0], nums[1]), ops[2][1](nums[2], nums[3]))
        result4 = ops[2][1](ops[0][1](nums[0], nums[1]), ops[1][1](nums[2], nums[3]))
        
        # Check if any result equals 24
        if any(abs(result - 24) < 1e-6 for result in [result1, result2, result3, result4]):
            return True
    except ZeroDivisionError:
        return False
    return False

def find_expression(nums):
    """Find an expression that evaluates to 24."""
    for num_perm in permutations(nums):
        for op_perm in permutations(ops, 3):
            if evaluate_expression(num_perm, op_perm):
                # Construct the expression string
                expr = f"(({num_perm[0]}{op_perm[0][0]}{num_perm[1]}){op_perm[1][0]}({num_perm[2]}{op_perm[2][0]}{num_perm[3]}))"
                return expr
    return None

# Define the numbers
numbers = [2, 2, 8, 12]

# Find and print the expression
expression = find_expression(numbers)
if expression:
    print(f"<<<{expression}>>>")
else:
    print("No solution found.")