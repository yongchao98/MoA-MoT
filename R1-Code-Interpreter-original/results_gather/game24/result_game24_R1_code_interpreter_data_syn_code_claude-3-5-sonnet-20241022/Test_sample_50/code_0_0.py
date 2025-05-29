from itertools import permutations, product
import re

def evaluate(expr):
    try:
        return eval(expr)
    except:
        return None

def find_24(numbers):
    # All possible operations
    ops = ['+', '-', '*', '/']
    
    # All possible number arrangements
    perms = list(permutations(numbers))
    
    # All possible operation combinations
    op_combs = list(product(ops, repeat=3))
    
    # Try all possible expressions with different parentheses patterns
    patterns = [
        "({0}{3}{1}){4}({2}{5}{6})",  # ((a?b)?(c?d))
        "({0}{3}({1}{4}{2})){5}{6}",  # (a?(b?c))?d
        "(({0}{3}{1}){4}{2}){5}{6}",  # ((a?b)?c)?d
        "{0}{3}({1}{4}({2}{5}{6}))",  # a?(b?(c?d))
    ]
    
    for nums in perms:
        for ops_set in op_combs:
            for pattern in patterns:
                expr = pattern.format(nums[0], nums[1], nums[2], ops_set[0], ops_set[1], ops_set[2], nums[3])
                result = evaluate(expr)
                if result is not None and abs(result - 24) < 1e-10:
                    return expr
    return None

numbers = [2, 12, 12, 12]
solution = find_24(numbers)
print(solution)