from itertools import permutations
from operator import add, sub, mul, truediv

def find_expression(numbers):
    # Define possible operations
    operations = [add, sub, mul, truediv]
    op_symbols = ['+', '-', '*', '/']

    def backtrack(nums, expr):
        if len(nums) == 1:
            if abs(nums[0] - 24) < 1e-6:  # Check if the result is 24
                return expr
            return None

        for i in range(len(nums)):
            for j in range(i + 1, len(nums)):
                for op, symbol in zip(operations, op_symbols):
                    # Try using nums[i] and nums[j] with each operation
                    try:
                        new_num = op(nums[i], nums[j])
                        new_expr = f"({expr[i]}{symbol}{expr[j]})"
                        new_nums = [nums[k] for k in range(len(nums)) if k != i and k != j] + [new_num]
                        new_exprs = [expr[k] for k in range(len(expr)) if k != i and k != j] + [new_expr]
                        result = backtrack(new_nums, new_exprs)
                        if result:
                            return result
                    except ZeroDivisionError:
                        continue
        return None

    # Try all permutations of numbers
    for perm in permutations(numbers):
        initial_expr = list(map(str, perm))
        result = backtrack(list(perm), initial_expr)
        if result:
            return f"<<<{result}>>>"

# Test the function with the given numbers
numbers = [4, 4, 11, 13]
print(find_expression(numbers))