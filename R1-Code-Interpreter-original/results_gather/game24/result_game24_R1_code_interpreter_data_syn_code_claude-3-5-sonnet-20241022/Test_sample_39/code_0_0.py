from itertools import permutations
import operator

def get_operations():
    return [
        (operator.add, '+'),
        (operator.sub, '-'),
        (operator.mul, '*'),
        (operator.truediv, '/')
    ]

def evaluate(nums, ops):
    # Different parentheses patterns
    patterns = [
        lambda n, o: o[0](o[1](o[2](n[0], n[1]), n[2]), n[3]),  # ((a?b)?c)?d
        lambda n, o: o[0](n[0], o[1](o[2](n[1], n[2]), n[3])),  # a?(b?(c?d))
        lambda n, o: o[0](o[1](n[0], n[1]), o[2](n[2], n[3])),  # (a?b)?(c?d)
    ]
    
    for pattern in patterns:
        try:
            result = pattern(nums, ops)
            if abs(result - 24) < 1e-10:  # Handle floating point precision
                return True
        except ZeroDivisionError:
            continue
    return False

def find_solution(numbers):
    ops = get_operations()
    
    # Try all number permutations
    for nums in permutations(numbers):
        # Try all operation combinations
        for op1 in ops:
            for op2 in ops:
                for op3 in ops:
                    if evaluate(nums, [op1[0], op2[0], op3[0]]):
                        # Found a solution, now reconstruct the expression
                        return f"(({nums[0]}{op1[1]}{nums[1]}){op2[1]}{nums[2]}){op3[1]}{nums[3]}"
    
    return "No solution found"

numbers = [2, 11, 11, 11]
print(find_solution(numbers))