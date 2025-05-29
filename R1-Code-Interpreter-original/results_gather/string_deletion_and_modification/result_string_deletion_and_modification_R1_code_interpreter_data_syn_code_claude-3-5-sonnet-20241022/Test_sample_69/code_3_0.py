def apply_rules(s):
    steps = [s]
    print(f"Starting string: {s}")
    
    while True:
        changed = False
        original = s
        
        # Rule 1: If starts with 'bb', remove second character
        if s.startswith('bb'):
            s = s[0] + s[2:]
            changed = True
            print(f"Rule 1: {s}")
        
        # Rule 2: If length > 15, remove middle character
        elif len(s) > 15:
            mid = len(s) // 2
            s = s[:mid] + s[mid+1:]
            changed = True
            print(f"Rule 2: {s}")
        
        # Rule 3: If ends with 'aa', replace with 'cc'
        elif s.endswith('aa'):
            s = s[:-2] + 'cc'
            changed = True
            print(f"Rule 3: {s}")
        
        # Rule 4: If starts with 'bc', delete first two and append 'aa'
        elif s.startswith('bc'):
            s = s[2:] + 'aa'
            changed = True
            print(f"Rule 4: {s}")
        
        # Rule 5: If even number of 'b's (and at least one), append 'ab'
        elif s.count('b') > 0 and s.count('b') % 2 == 0:
            s = s + 'ab'
            changed = True
            print(f"Rule 5: {s}")
        
        # Rule 6: If ends with 'bb', delete last two characters
        elif s.endswith('bb'):
            s = s[:-2]
            changed = True
            print(f"Rule 6: {s}")
        
        if changed:
            steps.append(s)
        else:
            break
    
    print("\nAll steps:")
    for i, step in enumerate(steps):
        print(f"Step {i}: {step}")
    print(f"\nFinal string: {s}")

s = "ababcaccbbbccbaabcbc"
apply_rules(s)