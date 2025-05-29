from itertools import permutations, product

def find_expression(numbers, target=24):
    # Helper function to evaluate expressions
    def evaluate(nums, expr):
        if len(nums) == 1:
            if abs(nums[0] - target) < 1e-6:  # Check if the result is close to the target
                print(f"<<<{expr}>>>")
                return True
            return False
        
        for i in range(len(nums)):
            for j in range(i + 1, len(nums)):
                # Try all operations
                for op in ['+', '-', '*', '/']:
                    if op == '/' and nums[j] == 0:
                        continue  # Avoid division by zero
                    
                    # Calculate new number
                    if op == '+':
                        new_num = nums[i] + nums[j]
                    elif op == '-':
                        new_num = nums[i] - nums[j]
                    elif op == '*':
                        new_num = nums[i] * nums[j]
                    elif op == '/':
                        new_num = nums[i] / nums[j]
                    
                    # Create new list of numbers
                    new_nums = [nums[k] for k in range(len(nums)) if k != i and k != j] + [new_num]
                    new_expr = f"({expr[i]}{op}{expr[j]})"
                    
                    # Recursively evaluate
                    if evaluate(new_nums, [new_expr] + [expr[k] for k in range(len(expr)) if k != i and k != j]):
                        return True
        return False

    # Try all permutations of numbers
    for perm in permutations(numbers):
        if evaluate(list(perm), list(map(str, perm))):
            break

# Example usage
numbers = [1, 5, 8, 9]
find_expression(numbers)