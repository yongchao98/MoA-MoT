def apply_operations(s):
    steps = [s]
    
    while True:
        original = s
        changed = False
        
        # Rule 5: contains 'bca' (special handling for 'bccba' pattern)
        i = 0
        while i < len(s):
            if i + 4 < len(s) and s[i] == 'b' and s[i+1] == 'c' and s[i+2] == 'c' and s[i+3] == 'b' and s[i+4] == 'a':
                # This handles the special case where 'bccba' contains 'bca'
                s = s[:i] + s[i+3:]  # Remove 'bcc'
                steps.append(f"{s} (Rule 5 - special case)")
                changed = True
                break
            elif i + 2 < len(s) and s[i] == 'b' and s[i+1] == 'c' and s[i+2] == 'a':
                s = s[:i] + s[i+3:]
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

    print("Initial string:", steps[0])
    print("\nSteps:")
    for i, step in enumerate(steps):
        print(f"Step {i}: {step}")
    print("\nFinal string:", s)

# Test with the string
s = "baccbabbcc"
apply_operations(s)