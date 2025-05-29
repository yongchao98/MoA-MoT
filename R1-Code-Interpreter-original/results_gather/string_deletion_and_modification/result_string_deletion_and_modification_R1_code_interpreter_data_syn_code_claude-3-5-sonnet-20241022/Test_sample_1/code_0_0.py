def operation1(s):
    # If ends with 'aa', replace with 'cc'
    if s.endswith('aa'):
        return s[:-2] + 'cc', True
    return s, False

def operation2(s):
    # If length > 15, remove middle character
    if len(s) > 15:
        mid = len(s) // 2
        return s[:mid] + s[mid+1:], True
    return s, False

def operation3(s):
    # If prefix is 'bc', delete first two chars and append 'aa'
    if s.startswith('bc'):
        return s[2:] + 'aa', True
    return s, False

def operation4(s):
    # If contains 'bca', delete first occurrence
    if 'bca' in s:
        idx = s.index('bca')
        return s[:idx] + s[idx+3:], True
    return s, False

def operation5(s):
    # If prefix is 'aa', remove first character
    if s.startswith('aa'):
        return s[1:], True
    return s, False

def operation6(s):
    # If starts with 'ac', replace first two chars with 'zz'
    if s.startswith('ac'):
        return 'zz' + s[2:], True
    return s, False

def process_string(s):
    operations = [operation1, operation2, operation3, operation4, operation5, operation6]
    changed = True
    print(f"Initial string: {s}")
    
    while changed:
        changed = False
        for op in operations:
            s, op_changed = op(s)
            if op_changed:
                print(f"After {op.__name__}: {s}")
                changed = True
                break
    
    return s

initial = "acaaaababcaccbbaaab"
final = process_string(initial)
print(f"\nFinal string: {final}")