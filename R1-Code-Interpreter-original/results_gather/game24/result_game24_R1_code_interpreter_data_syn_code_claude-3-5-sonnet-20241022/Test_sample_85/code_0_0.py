from itertools import permutations, product
import operator

def evaluate(nums, ops, paren):
    if paren == 0:  # ((a#b)#c)#d
        try:
            val1 = ops[0](nums[0], nums[1])
            val2 = ops[1](val1, nums[2])
            return ops[2](val2, nums[3])
        except:
            return None
    elif paren == 1:  # (a#b)#(c#d)
        try:
            val1 = ops[0](nums[0], nums[1])
            val2 = ops[2](nums[2], nums[3])
            return ops[1](val1, val2)
        except:
            return None
    elif paren == 2:  # (a#(b#c))#d
        try:
            val1 = ops[1](nums[1], nums[2])
            val2 = ops[0](nums[0], val1)
            return ops[2](val2, nums[3])
        except:
            return None

def find_24(numbers):
    operators = {'+': operator.add, '-': operator.sub, '*': operator.mul, '/': operator.truediv}
    ops_list = list(operators.values())
    
    for nums in permutations(numbers):
        for ops in product(ops_list, repeat=3):
            for paren in range(3):
                result = evaluate(nums, ops, paren)
                if result is not None and abs(result - 24) < 1e-10:
                    # Convert back to string representation
                    ops_str = []
                    for op in ops:
                        if op == operator.add: ops_str.append('+')
                        elif op == operator.sub: ops_str.append('-')
                        elif op == operator.mul: ops_str.append('*')
                        else: ops_str.append('/')
                    
                    if paren == 0:
                        return f"(({nums[0]}{ops_str[0]}{nums[1]}){ops_str[1]}{nums[2]}){ops_str[2]}{nums[3]}"
                    elif paren == 1:
                        return f"({nums[0]}{ops_str[0]}{nums[1]}){ops_str[1]}({nums[2]}{ops_str[2]}{nums[3]})"
                    else:
                        return f"({nums[0]}{ops_str[0]}({nums[1]}{ops_str[1]}{nums[2]})){ops_str[2]}{nums[3]}"
    
    return "No solution found"

print(find_24([1, 7, 10, 12]))