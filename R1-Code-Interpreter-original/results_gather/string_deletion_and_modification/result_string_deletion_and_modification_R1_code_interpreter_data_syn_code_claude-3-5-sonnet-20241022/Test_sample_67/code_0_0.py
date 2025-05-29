def operation1(s):
    # If starts with 'bb', remove second character
    if s.startswith('bb'):
        return s[0] + s[2:], True
    return s, False

def operation2(s):
    # If contains 'acb', replace first occurrence with 'bca'
    if 'acb' in s:
        idx = s.index('acb')
        return s[:idx] + 'bca' + s[idx+3:], True
    return s, False

def operation3(s):
    # If ends with 'aa', replace with 'cc'
    if s.endswith('aa'):
        return s[:-2] + 'cc', True
    return s, False

def operation4(s):
    # If length > 15, remove middle character
    if len(s) > 15:
        mid = len(s) // 2
        return s[:mid] + s[mid+1:], True
    return s, False

def operation5(s):
    # If ends with 'cc', replace with 'b' and prepend 'a'
    if s.endswith('cc'):
        return 'a' + s[:-2] + 'b', True
    return s, False

def operation6(s):
    # If starts with 'aa', remove first character
    if s.startswith('aa'):
        return s[1:], True
    return s, False

def process_string(s):
    operations = [operation1, operation2, operation3, operation4, operation5, operation6]
    changes_made = True
    steps = []
    steps.append(s)
    
    while changes_made:
        changes_made = False
        for op in operations:
            s, changed = op(s)
            if changed:
                steps.append(s)
                changes_made = True
                break
    
    print("Steps:")
    for i, step in enumerate(steps):
        print(f"Step {i}: {step}")
    return s

initial = "baacbcaccabcacbcbaac"
final = process_string(initial)
print(f"\nFinal result: {final}")