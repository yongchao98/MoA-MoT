def rule1(s):
    # Count b's and if even and at least one, append 'ab'
    if s.count('b') >= 1 and s.count('b') % 2 == 0:
        return s + 'ab', True
    return s, False

def rule2(s):
    # If prefix is 'cb', replace with 'aa' and delete last char
    if s.startswith('cb'):
        return 'aa' + s[2:-1], True
    return s, False

def rule3(s):
    # If contains 'bca', delete first occurrence
    if 'bca' in s:
        idx = s.index('bca')
        return s[:idx] + s[idx+3:], True
    return s, False

def rule4(s):
    # If ends with 'ca', remove last char
    if s.endswith('ca'):
        return s[:-1], True
    return s, False

def rule5(s):
    # If prefix is 'bc', delete first two chars and append 'aa'
    if s.startswith('bc'):
        return s[2:] + 'aa', True
    return s, False

def rule6(s):
    # If contains 'abc', replace first occurrence with 'cab'
    if 'abc' in s:
        idx = s.index('abc')
        return s[:idx] + 'cab' + s[idx+3:], True
    return s, False

def process_string(s):
    rules = [rule1, rule2, rule3, rule4, rule5, rule6]
    changed = True
    steps = []
    steps.append(s)
    
    while changed:
        changed = False
        for rule in rules:
            new_s, rule_applied = rule(s)
            if rule_applied:
                s = new_s
                steps.append(s)
                changed = True
                break
    
    print("\n".join(steps))
    return s

initial = "abbabbcaabcab"
final = process_string(initial)
print("\nFinal string:", final)