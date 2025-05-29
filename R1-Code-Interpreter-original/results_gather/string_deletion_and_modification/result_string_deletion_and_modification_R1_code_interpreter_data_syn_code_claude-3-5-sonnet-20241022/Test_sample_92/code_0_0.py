def rule1(s):
    # If 'ca' not at start, remove first occurrence after first char
    if len(s) > 2:
        pos = s[1:].find('ca')
        if pos != -1:
            return s[:pos+1] + s[pos+3:], True
    return s, False

def rule2(s):
    # Replace first 'acb' with 'bca'
    if 'acb' in s:
        pos = s.find('acb')
        return s[:pos] + 'bca' + s[pos+3:], True
    return s, False

def rule3(s):
    # If ends with 'aa', replace with 'cc'
    if s.endswith('aa'):
        return s[:-2] + 'cc', True
    return s, False

def rule4(s):
    # If even number of b's (and at least one), append 'ab'
    count_b = s.count('b')
    if count_b > 0 and count_b % 2 == 0:
        return s + 'ab', True
    return s, False

def rule5(s):
    # If starts with 'ca', replace with 'bb' and append 'c'
    if s.startswith('ca'):
        return 'bb' + s[2:] + 'c', True
    return s, False

def rule6(s):
    # If starts with 'cc', remove first two chars
    if s.startswith('cc'):
        return s[2:], True
    return s, False

def apply_rules(s):
    print(f"Starting with: {s}")
    iteration = 1
    while True:
        original = s
        changed = False
        
        # Try each rule in order
        for rule_num, rule in enumerate([rule1, rule2, rule3, rule4, rule5, rule6], 1):
            s, was_changed = rule(s)
            if was_changed:
                print(f"After rule {rule_num} (iteration {iteration}): {s}")
                changed = True
                break
        
        if not changed:
            break
        iteration += 1
    
    return s

initial = "acbbcabcccbcc"
final = apply_rules(initial)
print(f"\nFinal result: {final}")