def operation1(s):
    # If starts with 'ac', replace with 'zz'
    if s.startswith('ac'):
        return 'zz' + s[2:], True
    return s, False

def operation2(s):
    # If ends with 'bb', delete last two chars
    if s.endswith('bb'):
        return s[:-2], True
    return s, False

def operation3(s):
    # If even number of b's (and at least one), append 'ab'
    b_count = s.count('b')
    if b_count > 0 and b_count % 2 == 0:
        return s + 'ab', True
    return s, False

def operation4(s):
    # If starts with 'cb', replace with 'aa' and delete last char
    if s.startswith('cb'):
        return 'aa' + s[2:][:-1], True
    return s, False

def operation5(s):
    # If contains 'abc', replace first occurrence with 'cab'
    if 'abc' in s:
        idx = s.index('abc')
        return s[:idx] + 'cab' + s[idx+3:], True
    return s, False

def operation6(s):
    # If ends with 'aa', replace with 'cc'
    if s.endswith('aa'):
        return s[:-2] + 'cc', True
    return s, False

def process_string(s):
    operations = [operation1, operation2, operation3, operation4, operation5, operation6]
    changed = True
    steps = []
    steps.append(s)
    
    while changed:
        changed = False
        for op in operations:
            new_s, was_changed = op(s)
            if was_changed:
                s = new_s
                steps.append(s)
                changed = True
                break
    
    print("Steps:")
    for i, step in enumerate(steps):
        print(f"Step {i}: {step}")
    return s

initial = "accabcabccbcccccabca"
final = process_string(initial)
print(f"\nFinal result: {final}")