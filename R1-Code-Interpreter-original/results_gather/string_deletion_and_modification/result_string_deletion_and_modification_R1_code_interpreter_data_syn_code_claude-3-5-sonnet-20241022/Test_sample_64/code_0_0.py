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
    steps = 1
    while True:
        print(f"Step {steps}: {s}")
        original = s
        
        # Try each rule in sequence
        s = rule1(s)
        if s != original:
            continue
            
        s = rule2(s)
        if s != original:
            continue
            
        s = rule3(s)
        if s != original:
            continue
            
        s = rule4(s)
        if s != original:
            continue
            
        s = rule5(s)
        if s != original:
            continue
            
        s = rule6(s)
        if s != original:
            continue
            
        # If we get here, no rules applied
        break
        
        steps += 1
    
    return s

# Initial string
s = "baccabcaabcaabbcbca"
final = apply_rules(s)
print(f"\nFinal result: {final}")