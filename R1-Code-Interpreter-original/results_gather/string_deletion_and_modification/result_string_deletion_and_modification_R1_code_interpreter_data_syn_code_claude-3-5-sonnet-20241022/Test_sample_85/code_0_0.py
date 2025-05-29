def rule1(s):
    # If suffix is 'ac', replace with 'cb'
    if s.endswith('ac'):
        return s[:-2] + 'cb', True
    return s, False

def rule2(s):
    # If prefix is 'bc', delete first two chars and append 'aa'
    if s.startswith('bc'):
        return s[2:] + 'aa', True
    return s, False

def rule3(s):
    # If contains 'abc', replace first occurrence with 'cab'
    if 'abc' in s:
        idx = s.index('abc')
        return s[:idx] + 'cab' + s[idx+3:], True
    return s, False

def rule4(s):
    # If suffix is 'cc', replace with 'b' and prepend 'a'
    if s.endswith('cc'):
        return 'a' + s[:-2] + 'b', True
    return s, False

def rule5(s):
    # If ends with 'ca', remove last character
    if s.endswith('ca'):
        return s[:-1], True
    return s, False

def rule6(s):
    # If contains 'ca' (not at start), remove first occurrence after first char
    if 'ca' in s[1:]:
        idx = s[1:].index('ca') + 1
        return s[:idx] + s[idx+2:], True
    return s, False

def process_string(s):
    rules = [rule1, rule2, rule3, rule4, rule5, rule6]
    steps = [s]
    
    while True:
        original = s
        for rule in rules:
            s, changed = rule(s)
            if changed:
                steps.append(s)
                break
        if original == s:  # No rules applied
            break
    
    print("Initial string:", steps[0])
    for i, step in enumerate(steps[1:], 1):
        print(f"Step {i}: {step}")
    print("\nFinal string:", s)

# Test with the given string
input_string = "accacabbccbbbacaabc"
process_string(input_string)