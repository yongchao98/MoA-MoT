import itertools
import operator

# Define the numbers and target
numbers = [3, 8, 10, 11]
target = 24

# Define possible operations and their symbols
operations = {
    operator.add: '+',
    operator.sub: '-',
    operator.mul: '*',
    operator.truediv: '/'
}

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

def find_solution():
    # Generate all permutations of numbers
    for num_perm in itertools.permutations(numbers):
        # Generate all combinations of operations
        for ops in itertools.product(operations.keys(), repeat=3):
            results = evaluate_expression(num_perm, ops)
            for i, result in enumerate(results):
                if result == target:
                    # Format the expression based on the grouping
                    if i == 0:
                        expr = f"({num_perm[0]} {operations[ops[0]]} ({num_perm[1]} {operations[ops[1]]} ({num_perm[2]} {operations[ops[2]]} {num_perm[3]})))"
                    elif i == 1:
                        expr = f"(({num_perm[0]} {operations[ops[0]]} {num_perm[1]}) {operations[ops[1]]} ({num_perm[2]} {operations[ops[2]]} {num_perm[3]}))"
                    elif i == 2:
                        expr = f"(({num_perm[0]} {operations[ops[0]]} ({num_perm[1]} {operations[ops[1]]} {num_perm[2]})) {operations[ops[2]]} {num_perm[3]})"
                    elif i == 3:
                        expr = f"({num_perm[0]} {operations[ops[0]]} (({num_perm[1]} {operations[ops[1]]} {num_perm[2]}) {operations[ops[2]]} {num_perm[3]}))"
                    elif i == 4:
                        expr = f"(({num_perm[0]} {operations[ops[0]]} {num_perm[1]}) {operations[ops[1]]} ({num_perm[2]} {operations[ops[2]]} {num_perm[3]}))"
                    return f"<<<{expr} = {target}>>>"
    return "No solution found."

# Find and print the solution
print(find_solution())