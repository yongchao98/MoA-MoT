def operation1(s):
    # If suffix is 'bb', delete last two characters
    if s.endswith('bb'):
        return s[:-2], True
    return s, False

def operation2(s):
    # If contains 'ca' (not at start), remove first occurrence after first char
    if len(s) > 2:
        pos = s[1:].find('ca')
        if pos != -1:
            pos += 1  # Adjust position to account for the slice
            return s[:pos] + s[pos+2:], True
    return s, False

def operation3(s):
    # If suffix is 'cc', replace with 'b' and prepend 'a'
    if s.endswith('cc'):
        return 'a' + s[:-2] + 'b', True
    return s, False

def operation4(s):
    # If suffix is 'ac', replace with 'cb'
    if s.endswith('ac'):
        return s[:-2] + 'cb', True
    return s, False

def operation5(s):
    # If ends with 'ca', remove last character
    if s.endswith('ca'):
        return s[:-1], True
    return s, False

def operation6(s):
    # If prefix is 'cb', replace with 'aa' and delete last character
    if s.startswith('cb'):
        return 'aa' + s[2:-1], True
    return s, False

def process_string(s):
    operations = [operation1, operation2, operation3, operation4, operation5, operation6]
    changed = True
    steps = []
    steps.append(s)
    
    while changed:
        changed = False
        for op in operations:
            s, was_changed = op(s)
            if was_changed:
                steps.append(s)
                changed = True
                break
    
    print("Steps:")
    for i, step in enumerate(steps):
        print(f"Step {i}: {step}")
    return s

initial = "bbaaacccaaccabbca"
result = process_string(initial)
print(f"\nFinal result: {result}")