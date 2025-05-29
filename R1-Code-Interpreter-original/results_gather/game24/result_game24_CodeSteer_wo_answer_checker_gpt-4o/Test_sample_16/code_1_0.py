from itertools import permutations, product

def find_equation_to_24(numbers):
    # Define possible operations
    operations = ['+', '-', '*', '/']

    # Helper function to evaluate an expression
    def evaluate_expression(expr):
        try:
            return eval(expr)
        except ZeroDivisionError:
            return None

    # Recursive function to find the equation
    def search(nums, expr):
        if len(nums) == 1:
            if abs(nums[0] - 24) < 1e-6:  # Check if the result is 24
                print(f"<<<{expr} = 24>>>")
                return True
            return False

        for i in range(len(nums)):
            for j in range(i + 1, len(nums)):
                for op in operations:
                    # Form new numbers list
                    new_nums = nums[:i] + nums[i+1:j] + nums[j+1:]
                    if op == '+':
                        new_num = nums[i] + nums[j]
                    elif op == '-':
                        new_num = nums[i] - nums[j]
                    elif op == '*':
                        new_num = nums[i] * nums[j]
                    elif op == '/':
                        if nums[j] == 0:
                            continue
                        new_num = nums[i] / nums[j]

                    # Form new expression
                    new_expr = f"({expr[i]} {op} {expr[j]})"
                    if search(new_nums + [new_num], expr[:i] + expr[i+1:j] + expr[j+1:] + [new_expr]):
                        return True
        return False

    # Try all permutations of numbers
    for perm in permutations(numbers):
        if search(list(perm), list(map(str, perm))):
            return

# Example usage
numbers = [6, 7, 12, 12]
find_equation_to_24(numbers)