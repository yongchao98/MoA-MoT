from itertools import permutations

def apply_operations(nums, exprs, target=24):
    # Base case: if there's only one number left, check if it equals the target
    if len(nums) == 1:
        if abs(nums[0] - target) < 1e-6:  # Allow for floating-point precision issues
            return exprs[0]
        else:
            return None

    # Try all pairs of numbers and all operations
    for i in range(len(nums)):
        for j in range(len(nums)):
            if i != j:
                # Create a new list without the two numbers at i and j
                new_nums = [nums[k] for k in range(len(nums)) if k != i and k != j]
                new_exprs = [exprs[k] for k in range(len(exprs)) if k != i and k != j]

                # Try all operations
                for op in ['+', '-', '*', '/']:
                    if op == '+':
                        result = nums[i] + nums[j]
                        expression = f"({exprs[i]}+{exprs[j]})"
                    elif op == '-':
                        result = nums[i] - nums[j]
                        expression = f"({exprs[i]}-{exprs[j]})"
                    elif op == '*':
                        result = nums[i] * nums[j]
                        expression = f"({exprs[i]}*{exprs[j]})"
                    elif op == '/':
                        if nums[j] == 0:  # Avoid division by zero
                            continue
                        result = nums[i] / nums[j]
                        expression = f"({exprs[i]}/{exprs[j]})"

                    # Add the result to the new list and recurse
                    new_nums_with_result = new_nums + [result]
                    new_exprs_with_result = new_exprs + [expression]
                    sub_expression = apply_operations(new_nums_with_result, new_exprs_with_result, target)
                    if sub_expression:
                        return sub_expression

    return None

def find_expression_to_24(numbers):
    for perm in permutations(numbers):
        # Initialize expressions as strings of the numbers
        exprs = [str(num) for num in perm]
        expression = apply_operations(list(perm), exprs)
        if expression:
            return f"<<<{expression}>>>"
    return "No solution found."

# Test with the given input
numbers = [3, 7, 9, 9]
print(find_expression_to_24(numbers))