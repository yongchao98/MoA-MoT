def operation1(s):
    return (s[1:], True) if s.startswith('aa') else (s, False)

def operation2(s):
    b_count = s.count('b')
    return (s + 'ab', True) if b_count > 0 and b_count % 2 == 0 else (s, False)

def operation3(s):
    return ('ca' + s[2:], True) if s.startswith('ab') else (s, False)

def operation4(s):
    if 'acb' in s:
        idx = s.index('acb')
        return (s[:idx] + 'bca' + s[idx+3:], True)
    return (s, False)

def operation5(s):
    if len(s) > 1:
        rest = s[1:]
        if 'ca' in rest:
            idx = rest.index('ca')
            return (s[:idx+1] + s[idx+3:], True)
    return (s, False)

def operation6(s):
    if s.startswith('cb'):
        return ('aa' + s[2:-1], True)
    return (s, False)

def solve(s):
    operations = [operation1, operation2, operation3, operation4, operation5, operation6]
    steps = []
    steps.append(s)
    
    while True:
        original = s
        for op in operations:
            s, changed = op(s)
            if changed:
                steps.append(s)
                break
        if original == s:  # No changes made in full cycle
            break
    
    print("Initial:", steps[0])
    for i, step in enumerate(steps[1:], 1):
        print(f"Step {i}:", step)
    print("\nFinal string:", s)

# Test with the given input
solve("bbbcbbaaabbbb")