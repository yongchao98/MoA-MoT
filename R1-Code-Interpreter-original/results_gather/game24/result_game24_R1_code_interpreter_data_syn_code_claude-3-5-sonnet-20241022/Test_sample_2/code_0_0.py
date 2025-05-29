from itertools import permutations
import operator

def get_value(ops, nums):
    # All possible ways to place parentheses for 3 operations
    ways = [
        lambda a,b,c,d,op1,op2,op3: op1(op2(op3(a,b),c),d),
        lambda a,b,c,d,op1,op2,op3: op1(op2(a,op3(b,c)),d),
        lambda a,b,c,d,op1,op2,op3: op1(a,op2(op3(b,c),d)),
        lambda a,b,c,d,op1,op2,op3: op1(a,op2(b,op3(c,d)))
    ]
    
    operators = {
        '+': operator.add,
        '-': operator.sub,
        '*': operator.mul,
        '/': operator.truediv
    }
    
    # Try all possible combinations
    for nums_perm in permutations(nums):
        for op1 in operators.values():
            for op2 in operators.values():
                for op3 in operators.values():
                    for way in ways:
                        try:
                            result = way(*nums_perm, op1, op2, op3)
                            if abs(result - 24) < 1e-10:  # Account for floating point precision
                                return True, nums_perm, [op1, op2, op3], way
                        except ZeroDivisionError:
                            continue
    return False, None, None, None

def get_operator_symbol(op):
    if op == operator.add: return '+'
    if op == operator.sub: return '-'
    if op == operator.mul: return '*'
    if op == operator.truediv: return '/'

def format_expression(nums, ops, way_index):
    if way_index == 0:
        return f"({ops[0]} ({ops[1]} ({ops[2]} {nums[0]} {nums[1]}) {nums[2]}) {nums[3]})"
    elif way_index == 1:
        return f"({ops[0]} ({ops[1]} {nums[0]} ({ops[2]} {nums[1]} {nums[2]})) {nums[3]})"
    elif way_index == 2:
        return f"({ops[0]} {nums[0]} ({ops[1]} ({ops[2]} {nums[1]} {nums[2]}) {nums[3]}))"
    else:
        return f"({ops[0]} {nums[0]} ({ops[1]} {nums[1]} ({ops[2]} {nums[2]} {nums[3]})))"

numbers = [1, 1, 2, 8]
found, nums, ops, way = get_value(numbers)

if found:
    way_index = [lambda a,b,c,d,op1,op2,op3: op1(op2(op3(a,b),c),d),
                 lambda a,b,c,d,op1,op2,op3: op1(op2(a,op3(b,c)),d),
                 lambda a,b,c,d,op1,op2,op3: op1(a,op2(op3(b,c),d)),
                 lambda a,b,c,d,op1,op2,op3: op1(a,op2(b,op3(c,d)))].index(way)
    
    op_symbols = [get_operator_symbol(op) for op in ops]
    expression = format_expression(nums, op_symbols, way_index)
    print(expression)