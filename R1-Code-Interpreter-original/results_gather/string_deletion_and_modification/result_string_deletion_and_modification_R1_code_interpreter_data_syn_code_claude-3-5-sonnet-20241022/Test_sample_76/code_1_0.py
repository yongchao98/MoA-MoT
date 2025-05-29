def operation1(s):
    return (s[:-2], True) if s.endswith('bb') else (s, False)

def operation2(s):
    if len(s) > 2:
        pos = s[1:].find('ca')
        if pos != -1:
            pos += 1
            return (s[:pos] + s[pos+2:], True)
    return (s, False)

def operation3(s):
    return ('a' + s[:-2] + 'b', True) if s.endswith('cc') else (s, False)

def operation4(s):
    return (s[:-2] + 'cb', True) if s.endswith('ac') else (s, False)

def operation5(s):
    return (s[:-1], True) if s.endswith('ca') else (s, False)

def operation6(s):
    return ('aa' + s[2:-1], True) if s.startswith('cb') else (s, False)

def process_string(s):
    operations = [operation1, operation2, operation3, operation4, operation5, operation6]
    steps = [s]
    
    while True:
        original = s
        for op in operations:
            new_s, changed = op(s)
            if changed:
                s = new_s
                steps.append(s)
                break
        if original == s:  # No changes were made in the entire cycle
            break
    
    print("Steps:")
    for i, step in enumerate(steps):
        print(f"Step {i}: {step}")
    return s

initial = "bbaaacccaaccabbca"
result = process_string(initial)
print(f"\nFinal result: {result}")