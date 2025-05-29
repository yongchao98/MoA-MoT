import itertools
import operator

# Define the numbers and target
numbers = [3, 6, 6, 7]
target = 24

# Define the operations
operations = [('+', operator.add), ('-', operator.sub), ('*', operator.mul), ('/', operator.truediv)]

def evaluate_expression(nums, ops):
    try:
        # Apply the operations in sequence
        result = ops[0][1](nums[0], nums[1])
        result = ops[1][1](result, nums[2])
        result = ops[2][1](result, nums[3])
        return result
    except ZeroDivisionError:
        return None

# Try all permutations of numbers and combinations of operations
for num_perm in itertools.permutations(numbers):
    for ops in itertools.product(operations, repeat=3):
        if evaluate_expression(num_perm, ops) == target:
            # Format the answer
            answer = f"(({num_perm[0]}{ops[0][0]}{num_perm[1]}){ops[1][0]}{num_perm[2]}){ops[2][0]}{num_perm[3]}"
            print(f"<<<{answer}>>>")
            break