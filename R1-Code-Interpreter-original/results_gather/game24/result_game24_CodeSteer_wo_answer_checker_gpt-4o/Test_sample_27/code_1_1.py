from itertools import permutations
from sympy import symbols, simplify

def find_expression_to_24(numbers):
    def helper(nums, exprs):
        if len(nums) == 1:
            if abs(nums[0] - 24) < 1e-6:  # Check if the result is 24
                return exprs[0]
            return None
        
        for i in range(len(nums)):
            for j in range(i + 1, len(nums)):
                for op in ['+', '-', '*', '/']:
                    if op == '/' and nums[j] == 0:
                        continue  # Avoid division by zero
                    
                    new_nums = nums[:i] + nums[i+1:j] + nums[j+1:]
                    new_exprs = exprs[:i] + exprs[i+1:j] + exprs[j+1:]
                    
                    if op == '+':
                        new_num = nums[i] + nums[j]
                        new_expr = f"({exprs[i]}+{exprs[j]})"
                    elif op == '-':
                        new_num = nums[i] - nums[j]
                        new_expr = f"({exprs[i]}-{exprs[j]})"
                    elif op == '*':
                        new_num = nums[i] * nums[j]
                        new_expr = f"({exprs[i]}*{exprs[j]})"
                    elif op == '/':
                        new_num = nums[i] / nums[j]
                        new_expr = f"({exprs[i]}/{exprs[j]})"
                    
                    result = helper(new_nums + [new_num], new_exprs + [new_expr])
                    if result:
                        return result
        return None

    for perm in permutations(numbers):
        exprs = [str(num) for num in perm]
        result = helper(list(perm), exprs)
        if result:
            return result

numbers = [11, 13, 13, 13]
expression = find_expression_to_24(numbers)
print(f"<<<{expression}>>>")