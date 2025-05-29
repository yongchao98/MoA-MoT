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
    operators = {
        '+': operator.add,
        '-': operator.sub,
        '*': operator.mul,
        '/': operator.truediv
    }
    
    for nums in permutations(numbers):
        for ops in product([operators[op] for op in '+-*/'], repeat=3):
            for paren in range(3):
                result = evaluate(nums, ops, paren)
                if result is not None and abs(result - 24) < 1e-10:
                    # Convert operators back to symbols
                    op_symbols = []
                    for op in ops:
                        if op == operator.add: op_symbols.append('+')
                        elif op == operator.sub: op_symbols.append('-')
                        elif op == operator.mul: op_symbols.append('*')
                        else: op_symbols.append('/')
                    
                    if paren == 0:
                        return f"(({nums[0]}{op_symbols[0]}{nums[1]}){op_symbols[1]}{nums[2]}){op_symbols[2]}{nums[3]}"
                    elif paren == 1:
                        return f"({nums[0]}{op_symbols[0]}{nums[1]}){op_symbols[1]}({nums[2]}{op_symbols[2]}{nums[3]})"
                    else:
                        return f"({nums[0]}{op_symbols[0]}({nums[1]}{op_symbols[1]}{nums[2]})){op_symbols[2]}{nums[3]}"

result = find_24([1, 5, 6, 9])
print(result)