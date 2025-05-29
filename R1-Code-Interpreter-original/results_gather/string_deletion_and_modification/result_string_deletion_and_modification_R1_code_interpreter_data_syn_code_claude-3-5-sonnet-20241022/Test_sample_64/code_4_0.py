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
        
        # Try each rule in order
        # Rule 2 (length > 15)
        if len(s) > 15:
            mid = len(s) // 2
            s = s[:mid] + s[mid+1:]
            print(f"Applied Rule 2: removed character '{original[mid]}' at position {mid}")
            step += 1
            continue
        
        # Rule 1 (ends with 'cc')
        new_s = rule1(s)
        if new_s != s:
            print(f"Applied Rule 1: ends with 'cc'")
            s = new_s
            step += 1
            continue
            
        # Rule 3 (contains 'acb')
        new_s = rule3(s)
        if new_s != s:
            print(f"Applied Rule 3: 'acb' -> 'bca'")
            s = new_s
            step += 1
            continue
            
        # Rule 4 (starts with 'cb')
        new_s = rule4(s)
        if new_s != s:
            print(f"Applied Rule 4: starts with 'cb'")
            s = new_s
            step += 1
            continue
            
        # Rule 5 (starts with 'bb')
        new_s = rule5(s)
        if new_s != s:
            print(f"Applied Rule 5: starts with 'bb'")
            s = new_s
            step += 1
            continue
            
        # Rule 6 (starts with 'aa')
        new_s = rule6(s)
        if new_s != s:
            print(f"Applied Rule 6: starts with 'aa'")
            s = new_s
            step += 1
            continue
        
        if s == original:  # No rules applied
            break
    
    return s

initial = "baccabcaabcaabbcbca"
result = apply_rules(initial)
print(f"\nFinal result: {result}")