from itertools import permutations
from operator import add, sub, mul, truediv

def find_expression(nums):
    ops = [add, sub, mul, truediv]
    op_symbols = ['+', '-', '*', '/']

    def helper(nums):
        if len(nums) == 1:
            if abs(nums[0] - 24) < 1e-6:  # Check if the result is 24
                return str(nums[0])
            return None

        for i in range(len(nums)):
            for j in range(len(nums)):
                if i != j:
                    for op, symbol in zip(ops, op_symbols):
                        try:
                            new_nums = [nums[k] for k in range(len(nums)) if k != i and k != j]
                            new_nums.append(op(nums[i], nums[j]))
                            expr = helper(new_nums)
                            if expr:
                                return f"({nums[i]} {symbol} {nums[j]})" + expr
                        except ZeroDivisionError:
                            continue
        return None

    for perm in permutations(nums):
        result = helper(list(perm))
        if result:
            return result

nums = [5, 6, 12, 13]
expression = find_expression(nums)
if expression:
    print(f"<<<{expression}>>>")
else:
    print("No solution found.")