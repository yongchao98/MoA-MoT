def apply_rule1(s):
    # If contains 'bca', delete first occurrence
    if 'bca' in s:
        i = s.index('bca')
        return s[:i] + s[i+3:], True
    return s, False

def apply_rule2(s):
    # If suffix is 'cc', replace with 'b' and prepend 'a'
    if s.endswith('cc'):
        return 'a' + s[:-2] + 'b', True
    return s, False

def apply_rule3(s):
    # If starts with 'ac', replace with 'zz'
    if s.startswith('ac'):
        return 'zz' + s[2:], True
    return s, False

def apply_rule4(s):
    # If prefix is 'ab', replace with 'ca'
    if s.startswith('ab'):
        return 'ca' + s[2:], True
    return s, False

def apply_rule5(s):
    # If contains 'ca' (not at start), remove first such occurrence
    if len(s) > 1:
        i = s[1:].find('ca')
        if i != -1:
            i += 1  # adjust index to account for the slice
            return s[:i] + s[i+2:], True
    return s, False

def apply_rule6(s):
    # If prefix is 'bc', delete first two chars and append 'aa'
    if s.startswith('bc'):
        return s[2:] + 'aa', True
    return s, False

def process_string(s):
    print(f"Starting with: {s}")
    iteration = 1
    while True:
        original = s
        for rule_num, rule in enumerate([apply_rule1, apply_rule2, apply_rule3, 
                                       apply_rule4, apply_rule5, apply_rule6], 1):
            s, changed = rule(s)
            if changed:
                print(f"Iteration {iteration}, Rule {rule_num}: {s}")
                break
        if s == original:  # No rules applied
            break
        iteration += 1
    return s

final = process_string("bacccaccbca")
print(f"\nFinal result: {final}")