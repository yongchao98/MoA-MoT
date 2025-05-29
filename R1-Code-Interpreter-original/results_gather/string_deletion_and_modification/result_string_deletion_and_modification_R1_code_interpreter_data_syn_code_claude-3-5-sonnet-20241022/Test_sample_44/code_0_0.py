def operation1(s):
    # If ends with 'ca', remove last character
    if s.endswith('ca'):
        return s[:-1], True
    return s, False

def operation2(s):
    # If length > 15, remove middle character
    if len(s) > 15:
        mid = len(s) // 2
        return s[:mid] + s[mid+1:], True
    return s, False

def operation3(s):
    # If suffix is 'ac', replace with 'cb'
    if s.endswith('ac'):
        return s[:-2] + 'cb', True
    return s, False

def operation4(s):
    # If starts with 'bb', remove second character
    if s.startswith('bb'):
        return s[0] + s[2:], True
    return s, False

def operation5(s):
    # If contains 'acb', replace first occurrence with 'bca'
    if 'acb' in s:
        pos = s.find('acb')
        return s[:pos] + 'bca' + s[pos+3:], True
    return s, False

def operation6(s):
    # If starts with 'cc', remove first two characters
    if s.startswith('cc'):
        return s[2:], True
    return s, False

def process_string(s):
    operations = [operation1, operation2, operation3, operation4, operation5, operation6]
    changes_made = True
    
    print(f"Starting string: {s}")
    
    while changes_made:
        changes_made = False
        for op in operations:
            s, changed = op(s)
            if changed:
                print(f"After {op.__name__}: {s}")
                changes_made = True
                break
    
    return s

initial = "bcbbabccbcacbabacabb"
final = process_string(initial)
print(f"\nFinal string: {final}")