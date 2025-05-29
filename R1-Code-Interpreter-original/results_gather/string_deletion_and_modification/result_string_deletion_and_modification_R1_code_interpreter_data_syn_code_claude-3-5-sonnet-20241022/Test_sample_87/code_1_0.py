def apply_operations(s):
    steps = [s]
    
    while True:
        print(f"\nCurrent string: {s}")
        original = s
        changed = False
        
        # Rule 5: contains 'bca'
        # Important: in "baccba", when we remove one 'c', we get "bacba", which contains "bca"
        i = 0
        while i < len(s)-4:  # Need at least 5 characters for "baccba"
            if (s[i] == 'b' and 
                s[i+1] == 'a' and 
                s[i+2] == 'c' and 
                s[i+3] == 'c' and 
                s[i+4] == 'b' and 
                i+5 < len(s) and 
                s[i+5] == 'a'):
                # Found the pattern "baccba"
                print(f"Found pattern 'baccba' at position {i}")
                s = s[:i] + s[i+6:]
                steps.append(f"{s} (Rule 5)")
                changed = True
                break
            i += 1
        if changed:
            continue

        # Rule 1: 'ca' not at start
        if len(s) > 2 and 'ca' in s[1:]:
            pos = s[1:].find('ca') + 1
            s = s[:pos] + s[pos+2:]
            steps.append(f"{s} (Rule 1)")
            continue

        # Rule 2: ends with 'ba'
        if s.endswith('ba'):
            s = s[:-2] + 'ab'
            steps.append(f"{s} (Rule 2)")
            continue

        # Rule 3: ends with 'ac'
        if s.endswith('ac'):
            s = s[:-2] + 'cb'
            steps.append(f"{s} (Rule 3)")
            continue

        # Rule 4: starts with 'bb'
        if s.startswith('bb'):
            s = s[0] + s[2:]
            steps.append(f"{s} (Rule 4)")
            continue

        # Rule 6: starts with 'ca'
        if s.startswith('ca'):
            s = 'bb' + s[2:] + 'c'
            steps.append(f"{s} (Rule 6)")
            continue

        if s == original:
            break

    print("\nFinal steps:")
    for i, step in enumerate(steps):
        print(f"Step {i}: {step}")
    print("\nFinal string:", s)

# Test with the string
s = "baccbabbcc"
apply_operations(s)