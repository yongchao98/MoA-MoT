def operation1(s):
    # If ends with 'ba', replace with 'ab'
    if s.endswith('ba'):
        return s[:-2] + 'ab', True
    return s, False

def operation2(s):
    # If contains 'acb', replace first occurrence with 'bca'
    if 'acb' in s:
        idx = s.index('acb')
        return s[:idx] + 'bca' + s[idx+3:], True
    return s, False

def operation3(s):
    # If prefix is 'cb', replace with 'aa' and delete last char
    if s.startswith('cb'):
        return 'aa' + s[2:-1], True
    return s, False

def operation4(s):
    # If contains 'abc', replace first occurrence with 'cab'
    if 'abc' in s:
        idx = s.index('abc')
        return s[:idx] + 'cab' + s[idx+3:], True
    return s, False

def operation5(s):
    # If length > 15, remove middle character
    if len(s) > 15:
        mid = len(s) // 2
        return s[:mid] + s[mid+1:], True
    return s, False

def operation6(s):
    # If ends with 'aa', replace with 'cc'
    if s.endswith('aa'):
        return s[:-2] + 'cc', True
    return s, False

def process_string(s):
    operations = [operation1, operation2, operation3, operation4, operation5, operation6]
    step = 0
    print(f"Step {step}: {s}")
    
    while True:
        any_change = False
        for op_num, operation in enumerate(operations, 1):
            new_s, changed = operation(s)
            if changed:
                step += 1
                print(f"Step {step}: {new_s} (Operation {op_num})")
                s = new_s
                any_change = True
                break
        if not any_change:
            break
    
    return s

initial = "abaabbacbcabbaaaaaa"
final = process_string(initial)
print("\nFinal string:", final)