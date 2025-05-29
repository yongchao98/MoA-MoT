def operation1(s):
    # If contains 'bca', delete first occurrence
    if 'bca' in s:
        i = s.index('bca')
        return s[:i] + s[i+3:], True
    return s, False

def operation2(s):
    # If suffix is 'bb', delete last two chars
    if s.endswith('bb'):
        return s[:-2], True
    return s, False

def operation3(s):
    # If ends with 'ba', replace with 'ab'
    if s.endswith('ba'):
        return s[:-2] + 'ab', True
    return s, False

def operation4(s):
    # If ends with 'ca', remove last char
    if s.endswith('ca'):
        return s[:-1], True
    return s, False

def operation5(s):
    # If contains 'acb', replace first occurrence with 'bca'
    if 'acb' in s:
        i = s.index('acb')
        return s[:i] + 'bca' + s[i+3:], True
    return s, False

def operation6(s):
    # If starts with 'ac', replace with 'zz'
    if s.startswith('ac'):
        return 'zz' + s[2:], True
    return s, False

def process_string(s):
    operations = [operation1, operation2, operation3, operation4, operation5, operation6]
    changed = True
    steps = []
    
    while changed:
        changed = False
        for op in operations:
            new_s, was_changed = op(s)
            if was_changed:
                steps.append(f"{s} -> {new_s}")
                s = new_s
                changed = True
                break
    
    print("Final string:", s)
    print("\nSteps:")
    for step in steps:
        print(step)

# Test the string
input_string = "acacabababcacbcc"
process_string(input_string)