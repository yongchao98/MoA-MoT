def rule1(s):
    # Replace first 'acb' with 'bca'
    if 'acb' in s:
        i = s.index('acb')
        return s[:i] + 'bca' + s[i+3:], True
    return s, False

def rule2(s):
    # If starts with 'ab', replace with 'ca'
    if s.startswith('ab'):
        return 'ca' + s[2:], True
    return s, False

def rule3(s):
    # If ends with 'cc', replace with 'b' and prepend 'a'
    if s.endswith('cc'):
        return 'a' + s[:-2] + 'b', True
    return s, False

def rule4(s):
    # Remove first 'ca' after first character
    if len(s) > 1 and 'ca' in s[1:]:
        i = s[1:].index('ca') + 1
        return s[:i] + s[i+2:], True
    return s, False

def rule5(s):
    # Replace first 'abc' with 'cab'
    if 'abc' in s:
        i = s.index('abc')
        return s[:i] + 'cab' + s[i+3:], True
    return s, False

def rule6(s):
    # If starts with 'bb', remove second character
    if s.startswith('bb'):
        return s[0] + s[2:], True
    return s, False

def process_string(s):
    print(f"Starting with: {s}")
    iteration = 1
    while True:
        original = s
        for rule_num, rule in enumerate([rule1, rule2, rule3, rule4, rule5, rule6], 1):
            s, changed = rule(s)
            if changed:
                print(f"After rule {rule_num}: {s}")
                break
        if s == original:
            break
        iteration += 1
    return s

initial = "cbcccbcabcc"
final = process_string(initial)
print(f"\nFinal result: {final}")