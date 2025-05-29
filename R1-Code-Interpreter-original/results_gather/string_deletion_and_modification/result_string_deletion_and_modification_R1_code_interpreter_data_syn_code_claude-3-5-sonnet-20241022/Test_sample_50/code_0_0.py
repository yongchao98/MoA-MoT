def rule1(s):
    if s.endswith('cc'):
        return 'a' + s[:-2] + 'b'
    return s

def rule2(s):
    if s.endswith('ba'):
        return s[:-2] + 'ab'
    return s

def rule3(s):
    b_count = s.count('b')
    if b_count > 0 and b_count % 2 == 0:
        return s + 'ab'
    return s

def rule4(s):
    if s.startswith('bb'):
        return s[0] + s[2:]
    return s

def rule5(s):
    if len(s) > 1:
        idx = s[1:].find('ca')
        if idx != -1:
            idx += 1  # Adjust index to account for full string
            return s[:idx] + s[idx+2:]
    return s

def rule6(s):
    if s.startswith('ac'):
        return 'zz' + s[2:]
    return s

def apply_rules(s):
    steps = []
    steps.append(s)
    while True:
        original = s
        
        # Try each rule in order
        s = rule1(s)
        if s != original:
            steps.append(f"Rule 1: {s}")
            continue
            
        s = rule2(s)
        if s != original:
            steps.append(f"Rule 2: {s}")
            continue
            
        s = rule3(s)
        if s != original:
            steps.append(f"Rule 3: {s}")
            continue
            
        s = rule4(s)
        if s != original:
            steps.append(f"Rule 4: {s}")
            continue
            
        s = rule5(s)
        if s != original:
            steps.append(f"Rule 5: {s}")
            continue
            
        s = rule6(s)
        if s != original:
            steps.append(f"Rule 6: {s}")
            continue
            
        break
    
    print("\n".join(steps))
    return s

initial = "ccaacababc"
final = apply_rules(initial)
print("\nFinal result:", final)