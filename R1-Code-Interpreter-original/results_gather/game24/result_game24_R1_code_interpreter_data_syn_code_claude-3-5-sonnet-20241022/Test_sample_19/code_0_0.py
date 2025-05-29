from itertools import permutations
import operator

def evaluate(nums, ops):
    # Try different parentheses patterns
    try:
        # Pattern 1: ((a op1 b) op2 c) op3 d
        val1 = ops[0](nums[0], nums[1])
        val2 = ops[1](val1, nums[2])
        result1 = ops[2](val2, nums[3])
        
        # Pattern 2: (a op1 b) op2 (c op3 d)
        val3 = ops[0](nums[0], nums[1])
        val4 = ops[2](nums[2], nums[3])
        result2 = ops[1](val3, val4)
        
        # Pattern 3: a op1 (b op2 (c op3 d))
        val5 = ops[2](nums[2], nums[3])
        val6 = ops[1](nums[1], val5)
        result3 = ops[0](nums[0], val6)
        
        results = [result1, result2, result3]
        for r in results:
            if abs(r - 24) < 0.0001:  # Handle floating point comparison
                return True
    except:
        return False
    return False

def find_24(numbers):
    operators = [operator.add, operator.sub, operator.mul, operator.truediv]
    
    # Try all possible number arrangements
    for nums in permutations(numbers):
        # Try all possible operator combinations
        for op1 in operators:
            for op2 in operators:
                for op3 in operators:
                    if evaluate(nums, [op1, op2, op3]):
                        # Convert operators to symbols for display
                        op_symbols = []
                        for op in [op1, op2, op3]:
                            if op == operator.add: op_symbols.append('+')
                            elif op == operator.sub: op_symbols.append('-')
                            elif op == operator.mul: op_symbols.append('*')
                            else: op_symbols.append('/')
                        
                        # Pattern 1: ((a op1 b) op2 c) op3 d
                        print(f"((({nums[0]} {op_symbols[0]} {nums[1]}) {op_symbols[1]} {nums[2]}) {op_symbols[2]} {nums[3]}) = 24")
                        return

numbers = [4, 4, 11, 13]
find_24(numbers)