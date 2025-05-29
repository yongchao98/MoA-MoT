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
        if len(s) > 15:  # Always prioritize Rule 2 when length > 15
            mid = len(s) // 2
            s = s[:mid] + s[mid+1:]
            print(f"Applied Rule 2: removed character '{original[mid]}' at position {mid}")
            step += 1
            continue
            
        # Check for Rule 3 (acb -> bca)
        if 'acb' in s:
            pos = s.find('acb')
            s = s[:pos] + 'bca' + s[pos+3:]
            print(f"Applied Rule 3: replaced 'acb' with 'bca' at position {pos}")
            step += 1
            continue
            
        # Check for Rule 4 (cb -> aa)
        if s.startswith('cb'):
            s = 'aa' + s[2:-1]
            print(f"Applied Rule 4: replaced starting 'cb' with 'aa'")
            step += 1
            continue
            
        # Check for Rule 5 (bb)
        if s.startswith('bb'):
            s = s[0] + s[2:]
            print(f"Applied Rule 5: removed second character")
            step += 1
            continue
            
        # Check for Rule 6 (aa)
        if s.startswith('aa'):
            s = s[1:]
            print(f"Applied Rule 6: removed first character")
            step += 1
            continue
            
        # Check for Rule 1 (cc)
        if s.endswith('cc'):
            s = 'a' + s[:-2] + 'b'
            print(f"Applied Rule 1: replaced ending 'cc' with 'b' and prepended 'a'")
            step += 1
            continue
        
        if s == original:  # No rules applied
            break
    
    return s

initial = "baccabcaabcaabbcbca"
result = apply_rules(initial)
print(f"\nFinal result: {result}")