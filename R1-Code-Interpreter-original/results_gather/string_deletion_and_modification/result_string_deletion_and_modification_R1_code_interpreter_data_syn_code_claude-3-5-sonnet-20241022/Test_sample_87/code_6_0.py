def apply_operations(s):
    steps = [s]
    
    def find_bca(s):
        # For each 'b', look ahead for 'a' and check if all characters in between are 'c'
        for i in range(len(s)):
            if s[i] == 'b':
                for j in range(i+1, len(s)):
                    if s[j] == 'a':
                        # Check if all characters between i and j are 'c'
                        middle = s[i+1:j]
                        if all(c == 'c' for c in middle):
                            return (i, j+1)  # Return the range to remove
        return None

    while True:
        print(f"\nCurrent string: {s}")
        original = s
        
        # Rule 5: contains 'bca'
        bca_range = find_bca(s)
        if bca_range:
            start, end = bca_range
            print(f"Found 'bca' pattern from index {start} to {end}: {s[start:end]}")
            s = s[:start] + s[end:]
            steps.append(f"{s} (Rule 5: removed {s[start:end]})")
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