def rule1(s):
    if s.endswith('cc'):
        return 'a' + s[:-2] + 'b'
    return s

def rule2(s):
    if len(s) > 15:
        mid = len(s) // 2
        return s[:mid] + s[mid+1:]
    return s

def rule3(s):
    if 'acb' in s:
        pos = s.find('acb')
        return s[:pos] + 'bca' + s[pos+3:]
    return s

def rule4(s):
    if s.startswith('cb'):
        return 'aa' + s[2:-1]
    return s

def rule5(s):
    if s.startswith('bb'):
        return s[0] + s[2:]
    return s

def rule6(s):
    if s.startswith('aa'):
        return s[1:]
    return s

def apply_rules(s):
    step = 1
    while True:
        print(f"Step {step}: {s} (length: {len(s)})")
        original = s
        
        # Try each rule
        for rule_num, rule_func in enumerate([rule1, rule2, rule3, rule4, rule5, rule6], 1):
            new_s = rule_func(s)
            if new_s != s:
                print(f"Applied Rule {rule_num}: {new_s}")
                s = new_s
                step += 1
                break
        
        if s == original:  # No rules applied
            break
    
    return s

initial = "baccabcaabcaabbcbca"
result = apply_rules(initial)
print(f"\nFinal result: {result}")