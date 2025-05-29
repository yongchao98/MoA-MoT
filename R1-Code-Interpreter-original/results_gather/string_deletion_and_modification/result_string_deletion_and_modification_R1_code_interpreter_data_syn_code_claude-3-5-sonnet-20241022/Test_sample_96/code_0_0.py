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
    # If contains 'ca' (not at start), remove first occurrence after first char
    if len(s) > 1 and 'ca' in s[1:]:
        idx = s.index('ca', 1)
        return s[:idx] + s[idx+2:], True
    return s, False

def operation4(s):
    # If contains even number of 'b's (and at least one), append 'ab'
    b_count = s.count('b')
    if b_count > 0 and b_count % 2 == 0:
        return s + 'ab', True
    return s, False

def operation5(s):
    # If prefix is 'ab', replace with 'ca'
    if s.startswith('ab'):
        return 'ca' + s[2:], True
    return s, False

def operation6(s):
    # If prefix is 'ca', replace with 'bb' and append 'c'
    if s.startswith('ca'):
        return 'bb' + s[2:] + 'c', True
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
                steps.append(s)
                break
        if original == s:
            break
    
    print(f"Final string: {s}")
    print("\nSteps:")
    for i, step in enumerate(steps):
        print(f"Step {i}: {step}")

# Test with the given string
input_string = "abacbbbcbbbacbcacba"
process_string(input_string)