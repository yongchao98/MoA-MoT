def find_equation_to_24(numbers):
    def helper(nums, equation):
        if len(nums) == 1:
            if abs(nums[0] - 24) < 1e-6:  # Check if the result is 24
                return equation
            return None
        
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
                    
                    # Form new list of numbers
                    new_nums = [nums[k] for k in range(len(nums)) if k != i and k != j]
                    new_nums.append(new_num)
                    
                    # Form new equation
                    new_equation = f"({equation[i]} {op} {equation[j]})"
                    
                    # Recursive call
                    result = helper(new_nums, [new_equation] + [equation[k] for k in range(len(equation)) if k != i and k != j])
                    if result:
                        return result
        return None

    # Initial call with numbers and their string representations
    initial_equation = [str(num) for num in numbers]
    return helper(numbers, initial_equation)

# Test the function with the given input
numbers = [1, 2, 3, 12]
equation = find_equation_to_24(numbers)
print(f"<<<{equation}>>>")