def rule1(s):
    if len(s) > 2:
        pos = s[1:].find('ca')
        if pos != -1:
            return s[:pos+1] + s[pos+3:], True
    return s, False

def rule2(s):
    if 'acb' in s:
        pos = s.find('acb')
        return s[:pos] + 'bca' + s[pos+3:], True
    return s, False

def rule3(s):
    if s.endswith('aa'):
        return s[:-2] + 'cc', True
    return s, False

def rule4(s):
    count_b = s.count('b')
    if count_b > 0 and count_b % 2 == 0:
        return s + 'ab', True
    return s, False

def rule5(s):
    if s.startswith('ca'):
        return 'bb' + s[2:] + 'c', True
    return s, False

def rule6(s):
    if s.startswith('cc'):
        return s[2:], True
    return s, False

def apply_rules(s):
    print(f"Starting with: {s}")
    iteration = 1
    while True:
        original = s
        for rule_num, rule in enumerate([rule1, rule2, rule3, rule4, rule5, rule6], 1):
            s, changed = rule(s)
            if changed:
                print(f"After rule {rule_num} (iteration {iteration}): {s}")
                iteration += 1
                break
        if original == s:
            break
    return s

initial = "acbbcabcccbcc"
final = apply_rules(initial)
print(f"\nFinal result: {final}")