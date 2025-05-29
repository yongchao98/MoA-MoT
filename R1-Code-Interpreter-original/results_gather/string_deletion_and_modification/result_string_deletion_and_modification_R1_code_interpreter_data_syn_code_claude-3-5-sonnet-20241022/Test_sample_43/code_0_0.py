def operation1(s):
    # If suffix is 'cc', replace with 'b' and prepend 'a'
    if s.endswith('cc'):
        return 'a' + s[:-2] + 'b'
    return s

def operation2(s):
    # If starts with 'bb', remove second character
    if s.startswith('bb'):
        return s[0] + s[2:]
    return s

def operation3(s):
    # If ends with 'aa', replace with 'cc'
    if s.endswith('aa'):
        return s[:-2] + 'cc'
    return s

def operation4(s):
    # If length > 15, remove middle character
    if len(s) > 15:
        mid = len(s) // 2
        return s[:mid] + s[mid+1:]
    return s

def operation5(s):
    # If suffix is 'ac', replace with 'cb'
    if s.endswith('ac'):
        return s[:-2] + 'cb'
    return s

def operation6(s):
    # If even number of 'b's (and at least one), append 'ab'
    count_b = s.count('b')
    if count_b > 0 and count_b % 2 == 0:
        return s + 'ab'
    return s

def apply_operations(s):
    operations = [operation1, operation2, operation3, operation4, operation5, operation6]
    steps = []
    steps.append(s)
    
    while True:
        changed = False
        for op in operations:
            new_s = op(s)
            if new_s != s:
                s = new_s
                steps.append(s)
                changed = True
                break
        if not changed:
            break
    
    print("Steps:")
    for i, step in enumerate(steps):
        print(f"Step {i}: {step}")
    return s

initial = "aabbabaacbaabba"
result = apply_operations(initial)
print(f"\nFinal result: {result}")