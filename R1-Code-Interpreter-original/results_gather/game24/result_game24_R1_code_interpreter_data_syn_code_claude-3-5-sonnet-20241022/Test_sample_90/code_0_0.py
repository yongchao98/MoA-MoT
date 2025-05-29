from itertools import permutations
import operator

def get_value(ops, nums):
    # All possible ways to parenthesize 3 operations
    ways = [
        lambda a,b,c,d: (a+b+c)+d,
        lambda a,b,c,d: (a+b)+(c+d),
        lambda a,b,c,d: a+(b+c+d),
        lambda a,b,c,d: ((a+b)+c)+d,
        lambda a,b,c,d: (a+(b+c))+d,
        lambda a,b,c,d: a+((b+c)+d),
        lambda a,b,c,d: a+(b+(c+d))
    ]
    
    operators = {
        '+': operator.add,
        '-': operator.sub,
        '*': operator.mul,
        '/': operator.truediv
    }
    
    results = []
    for way in ways:
        try:
            val = way(*[operators[ops[i]](nums[i], nums[i+1]) for i in range(3)])
            if abs(val - 24) < 0.0001:
                return True
        except:
            continue
    return False

def find_24(numbers):
    ops = ['+', '-', '*', '/']
    for nums in permutations(numbers):
        for op1 in ops:
            for op2 in ops:
                for op3 in ops:
                    if get_value([op1, op2, op3], nums):
                        # Try some common patterns
                        if op1 == '*' and op2 == '+' and op3 == '/':
                            expr = f"(({nums[0]}*{nums[1]})+{nums[2]})/{nums[3]}"
                            if abs(eval(expr) - 24) < 0.0001:
                                print(f"{expr} = 24")
                                return
                        if op1 == '*' and op2 == '-' and op3 == '/':
                            expr = f"(({nums[0]}*{nums[1]})-{nums[2]})/{nums[3]}"
                            if abs(eval(expr) - 24) < 0.0001:
                                print(f"{expr} = 24")
                                return
                        if op1 == '*' and op2 == '*' and op3 == '/':
                            expr = f"({nums[0]}*{nums[1]}*{nums[2]})/{nums[3]}"
                            if abs(eval(expr) - 24) < 0.0001:
                                print(f"{expr} = 24")
                                return

numbers = [4, 6, 6, 8]
find_24(numbers)