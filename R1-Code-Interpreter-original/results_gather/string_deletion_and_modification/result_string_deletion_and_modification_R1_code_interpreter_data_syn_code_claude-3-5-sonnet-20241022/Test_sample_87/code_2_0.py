def apply_operations(s):
    steps = [s]
    
    while True:
        changed = False
        current = s
        
        # Rule 1: 'ca' not at start
        if len(s) > 2 and 'ca' in s[1:]:
            pos = s[1:].find('ca') + 1
            s = s[:pos] + s[pos+2:]
            steps.append(f"{s} (Rule 1)")
            changed = True
        
        # Rule 2: ends with 'ba'
        elif s.endswith('ba'):
            s = s[:-2] + 'ab'
            steps.append(f"{s} (Rule 2)")
            changed = True
            
        # Rule 3: ends with 'ac'
        elif s.endswith('ac'):
            s = s[:-2] + 'cb'
            steps.append(f"{s} (Rule 3)")
            changed = True
            
        # Rule 4: starts with 'bb'
        elif s.startswith('bb'):
            s = s[0] + s[2:]
            steps.append(f"{s} (Rule 4)")
            changed = True
            
        # Rule 5: contains 'bca'
        elif 'bca' in s:
            pos = s.find('bca')
            s = s[:pos] + s[pos+3:]
            steps.append(f"{s} (Rule 5)")
            changed = True
            
        # Rule 6: starts with 'ca'
        elif s.startswith('ca'):
            s = 'bb' + s[2:] + 'c'
            steps.append(f"{s} (Rule 6)")
            changed = True
            
        if not changed:
            break
    
    print("Initial string:", steps[0])
    print("\nSteps:")
    for i, step in enumerate(steps):
        print(f"Step {i}: {step}")
    print("\nFinal string:", s)

# Test with the given string
s = "baccbabbcc"
apply_operations(s)