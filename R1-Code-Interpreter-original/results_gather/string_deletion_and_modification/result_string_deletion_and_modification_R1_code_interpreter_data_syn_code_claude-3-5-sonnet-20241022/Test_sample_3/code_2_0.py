def apply_rule1(s):
    return s[:-2] + 'cb' if s.endswith('ac') else s

def apply_rule2(s):
    if len(s) > 15:
        mid = len(s) // 2
        return s[:mid] + s[mid+1:]
    return s

def apply_rule3(s):
    return s[:-1] if s.endswith('ca') else s

def apply_rule4(s):
    b_count = s.count('b')
    return s + 'ab' if b_count > 0 and b_count % 2 == 0 else s

def apply_rule5(s):
    if 'bca' in s:
        idx = s.find('bca')
        return s[:idx] + s[idx+3:]
    return s

def apply_rule6(s):
    return 'zz' + s[2:] if s.startswith('ac') else s

def solve(s):
    current = s
    while True:
        changed = False
        old = current
        
        # Try each rule in sequence
        rules = [apply_rule1, apply_rule2, apply_rule3, apply_rule4, apply_rule5, apply_rule6]
        for rule in rules:
            new = rule(current)
            if new != current:
                current = new
                changed = True
                
        if not changed:
            break
            
    return current

initial = "cbaccccccbcabcbcabcb"
result = solve(initial)
print(f"Final result: {result}")