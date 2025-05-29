from itertools import permutations

def apply_operations(nums, target):
    # Helper function to recursively apply operations
    def helper(nums, expr):
        if len(nums) == 1:
            if abs(nums[0] - target) < 1e-6:  # Check if the result is close to the target
                return expr
            return None
        
        for i in range(len(nums) - 1):
            for op in ['+', '-', '*', '/']:
                if op == '/' and nums[i+1] == 0:
                    continue  # Avoid division by zero
                
                new_nums = nums[:i] + [eval(f"{nums[i]}{op}{nums[i+1]}")] + nums[i+2:]
                new_expr = f"({expr[i]}{op}{expr[i+1]})"
                result = helper(new_nums, expr[:i] + [new_expr] + expr[i+2:])
                if result:
                    return result
        return None

    # Try all permutations of the numbers
    for perm in permutations(nums):
        expr = list(map(str, perm))
        result = helper(list(perm), expr)
        if result:
            return result
    return "No solution found"

# Input numbers
numbers = [3, 6, 9, 9]
target = 24

# Find and print the solution
solution = apply_operations(numbers, target)
print(f"<<<{solution} = {target}>>>")