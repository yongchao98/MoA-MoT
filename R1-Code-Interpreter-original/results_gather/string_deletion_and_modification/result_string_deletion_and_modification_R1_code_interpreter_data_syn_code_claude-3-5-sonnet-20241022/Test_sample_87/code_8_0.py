def apply_operations(s):
    steps = [s]
    
    while True:
        print(f"\nCurrent string: {s}")
        original = s
        changed = False
        
        # Rule 5: contains 'bca'
        # Look for 'b' followed by any number of 'c's and then 'a'
        i = 0
        while i < len(s):
            if s[i] == 'b':
                j = i + 1
                # Count consecutive 'c's
                while j < len(s) and s[j] == 'c':
                    j += 1
                # Check if next character is 'a'
                if j < len(s) and s[j] == 'a':
                    # We found a pattern: b + (some c's) + a
                    if j > i + 1:  # Make sure there's at least one 'c'
                        print(f"Found pattern: {s[i:j+1]}")
                        s = s[:i] + s[j+1:]
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