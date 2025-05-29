def operation1(s):
    # If starts with 'cc', remove first two chars
    if s.startswith('cc'):
        return s[2:], True
    return s, False

def operation2(s):
    # If contains 'acb', replace first occurrence with 'bca'
    if 'acb' in s:
        pos = s.find('acb')
        return s[:pos] + 'bca' + s[pos+3:], True
    return s, False

def operation3(s):
    # If prefix is 'aa', remove first char
    if s.startswith('aa'):
        return s[1:], True
    return s, False

def operation4(s):
    # If ends with 'ca', remove last char
    if s.endswith('ca'):
        return s[:-1], True
    return s, False

def operation5(s):
    # If suffix is 'cc', replace with 'b' and prepend 'a'
    if s.endswith('cc'):
        return 'a' + s[:-2] + 'b', True
    return s, False

def operation6(s):
    # If even number of 'b's (and at least one), append 'ab'
    b_count = s.count('b')
    if b_count > 0 and b_count % 2 == 0:
        return s + 'ab', True
    return s, False

def solve(s):
    operations = [operation1, operation2, operation3, operation4, operation5, operation6]
    changed = True
    steps = []
    
    while changed:
        changed = False
        for op in operations:
            s, was_changed = op(s)
            if was_changed:
                steps.append(s)
                changed = True
                break
    
    print(f"Final string: {s}")
    print("\nSteps:")
    for i, step in enumerate(steps, 1):
        print(f"Step {i}: {step}")

# Test with the given string
initial = "aaabcbbccbcaacc"
solve(initial)