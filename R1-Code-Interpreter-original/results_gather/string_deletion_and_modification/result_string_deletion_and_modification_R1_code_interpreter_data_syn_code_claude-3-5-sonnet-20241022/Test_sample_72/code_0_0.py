def operation1(s):
    # Count b's and append 'ab' if even number of b's
    if s.count('b') >= 1 and s.count('b') % 2 == 0:
        return s + 'ab', True
    return s, False

def operation2(s):
    # Find 'ca' after first character and remove it
    if len(s) > 2:
        idx = s[1:].find('ca')
        if idx != -1:
            idx += 1  # Adjust index to account for the slice
            return s[:idx] + s[idx+2:], True
    return s, False

def operation3(s):
    # Remove first two chars if starts with 'cc'
    if s.startswith('cc'):
        return s[2:], True
    return s, False

def operation4(s):
    # Replace 'ac' with 'zz' at start
    if s.startswith('ac'):
        return 'zz' + s[2:], True
    return s, False

def operation5(s):
    # Remove first char if starts with 'aa'
    if s.startswith('aa'):
        return s[1:], True
    return s, False

def operation6(s):
    # Replace 'cb' with 'aa' at start and remove last char
    if s.startswith('cb'):
        return 'aa' + s[2:-1], True
    return s, False

def process_string(s):
    operations = [operation1, operation2, operation3, operation4, operation5, operation6]
    steps = []
    steps.append(s)
    
    while True:
        original = s
        for op in operations:
            s, changed = op(s)
            if changed:
                steps.append(f"{s} (after {op.__name__})")
                break
        if original == s:
            break
    
    print("\n".join(steps))
    return s

initial = "cbaacbcabaacb"
final = process_string(initial)
print(f"\nFinal result: {final}")