def apply_operations(s):
    steps = [s]
    
    def find_bca(s):
        # Check for 'bca' including when there are extra characters between b and ca
        for i in range(len(s)-1):
            if s[i] == 'b':
                # Look for 'ca' after 'b', allowing for characters in between
                for j in range(i+1, len(s)-1):
                    if s[j:j+2] == 'ca':
                        # Check if this forms a valid 'bca' pattern
                        between = s[i+1:j]
                        if all(c == 'c' for c in between):  # All characters between are 'c'
                            return i, j+2
        return -1, -1

    while True:
        original = s
        changed = False
        
        # Print current state and analysis
        print(f"Analyzing string: {s}")
        
        # Rule 5: contains 'bca' (including when separated by 'c's)
        start, end = find_bca(s)
        if start != -1:
            s = s[:start] + s[end:]
            steps.append(f"{s} (Rule 5)")
            print(f"Found bca pattern from index {start} to {end}")
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

    print("\nSteps:")
    for i, step in enumerate(steps):
        print(f"Step {i}: {step}")
    print("\nFinal string:", s)

# Test with the string
s = "baccbabbcc"
apply_operations(s)