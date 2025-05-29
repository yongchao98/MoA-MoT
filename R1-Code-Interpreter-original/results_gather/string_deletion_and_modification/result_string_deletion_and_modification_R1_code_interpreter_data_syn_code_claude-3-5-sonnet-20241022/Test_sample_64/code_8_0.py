def apply_rules(s):
    step = 1
    while True:
        print(f"\nStep {step}: {s}")
        
        # Store original string to check if any changes were made
        original = s
        
        # Check each character position for 'acb' pattern
        for i in range(len(s)-2):
            if s[i:i+3] == 'acb':
                print(f"Found 'acb' at position {i}: {s[:i]}|{s[i:i+3]}|{s[i+3:]}")
                s = s[:i] + 'bca' + s[i+3:]
                print(f"After Rule 3: {s}")
                break
                
        # If no 'acb' found, check other rules
        if s == original:
            # Rule 5: starts with 'bb'
            if s.startswith('bb'):
                print(f"Found 'bb' at start: {s[:2]}|{s[2:]}")
                s = s[0] + s[2:]
                print(f"After Rule 5: {s}")
            # Rule 4: starts with 'cb'
            elif s.startswith('cb'):
                print(f"Found 'cb' at start: {s[:2]}|{s[2:]}")
                s = 'aa' + s[2:-1]
                print(f"After Rule 4: {s}")
            # Rule 6: starts with 'aa'
            elif s.startswith('aa'):
                print(f"Found 'aa' at start: {s[:2]}|{s[2:]}")
                s = s[1:]
                print(f"After Rule 6: {s}")
            # Rule 1: ends with 'cc'
            elif s.endswith('cc'):
                print(f"Found 'cc' at end: {s[:-2]}|{s[-2:]}")
                s = 'a' + s[:-2] + 'b'
                print(f"After Rule 1: {s}")
        
        # If no changes were made, we're done
        if s == original:
            break
            
        step += 1
    
    return s

s = "baccabcaabbcbca"
result = apply_rules(s)
print(f"\nFinal result: {result}")