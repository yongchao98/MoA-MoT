def rule1(s):
    # Count b's and if even and at least one, append 'ab'
    b_count = s.count('b')
    if b_count > 0 and b_count % 2 == 0:
        return s + 'ab', True
    return s, False

def rule2(s):
    # If prefix is 'ab', replace with 'ca'
    if s.startswith('ab'):
        return 'ca' + s[2:], True
    return s, False

def rule3(s):
    # If suffix is 'ac', replace with 'cb'
    if s.endswith('ac'):
        return s[:-2] + 'cb', True
    return s, False

def rule4(s):
    # If starts with 'bb', remove second character
    if s.startswith('bb'):
        return s[0] + s[2:], True
    return s, False

def rule5(s):
    # If suffix is 'cc', replace with 'b' and prepend 'a'
    if s.endswith('cc'):
        return 'a' + s[:-2] + 'b', True
    return s, False

def rule6(s):
    # If length > 15, remove middle character
    if len(s) > 15:
        mid = len(s) // 2
        return s[:mid] + s[mid+1:], True
    return s, False

def apply_rules(s):
    rules = [rule1, rule2, rule3, rule4, rule5, rule6]
    steps = []
    steps.append(s)
    
    while True:
        changed = False
        for rule in rules:
            new_s, rule_changed = rule(s)
            if rule_changed:
                s = new_s
                steps.append(s)
                changed = True
                break
        if not changed:
            break
    
    print("Steps:")
    for i, step in enumerate(steps):
        print(f"Step {i}: {step}")
    return s

initial = "aabbbcbcaaaabaaab"
final = apply_rules(initial)
print(f"\nFinal result: {final}")