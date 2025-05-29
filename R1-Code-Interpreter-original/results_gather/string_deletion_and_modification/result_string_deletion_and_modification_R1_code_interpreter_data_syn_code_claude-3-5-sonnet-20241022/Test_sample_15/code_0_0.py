def rule1(s):
    if s.startswith('cb'):
        return 'aa' + s[2:-1]
    return s

def rule2(s):
    if len(s) > 15:
        mid = len(s) // 2
        return s[:mid] + s[mid+1:]
    return s

def rule3(s):
    if 'abc' in s:
        pos = s.find('abc')
        return s[:pos] + 'cab' + s[pos+3:]
    return s

def rule4(s):
    if s.endswith('bb'):
        return s[:-2]
    return s

def rule5(s):
    if s.endswith('aa'):
        return s[:-2] + 'cc'
    return s

def rule6(s):
    if s.endswith('ca'):
        return s[:-1]
    return s

def apply_rules(s):
    rules = [rule1, rule2, rule3, rule4, rule5, rule6]
    iteration = 1
    while True:
        original = s
        for i, rule in enumerate(rules, 1):
            new_s = rule(s)
            if new_s != s:
                print(f"Iteration {iteration}, Rule {i}: {s} -> {new_s}")
                s = new_s
                break
        if original == s:
            break
        iteration += 1
    return s

initial = "bccccbcccababcccacb"
final = apply_rules(initial)
print(f"\nFinal result: {final}")