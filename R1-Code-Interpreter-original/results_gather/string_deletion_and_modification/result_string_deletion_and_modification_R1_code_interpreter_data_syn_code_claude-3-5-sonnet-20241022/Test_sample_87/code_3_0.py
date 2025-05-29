def apply_operations(s):
    def find_all_occurrences(s, sub):
        start = 0
        while True:
            start = s.find(sub, start)
            if start == -1:
                break
            yield start
            start += 1

    steps = [s]
    while True:
        original = s
        changed = False
        
        # Debug prints to see what patterns we're finding
        print(f"Current string: {s}")
        print(f"Contains 'bca'? {list(find_all_occurrences(s, 'bca'))}")
        print(f"Contains 'ca' after first char? {list(find_all_occurrences(s[1:], 'ca'))}")
        print(f"Ends with 'ba'? {s.endswith('ba')}")
        print(f"Ends with 'ac'? {s.endswith('ac')}")
        print(f"Starts with 'bb'? {s.startswith('bb')}")
        print(f"Starts with 'ca'? {s.startswith('ca')}")
        print("---")

        # Rule 5: contains 'bca' (checking this first as it's most relevant for current string)
        if 'bca' in s:
            pos = s.find('bca')
            s = s[:pos] + s[pos+3:]
            steps.append(f"{s} (Rule 5)")
            changed = True
            continue

        # Rule 1: 'ca' not at start
        if len(s) > 2 and 'ca' in s[1:]:
            pos = s[1:].find('ca') + 1
            s = s[:pos] + s[pos+2:]
            steps.append(f"{s} (Rule 1)")
            changed = True
            continue

        # Rule 2: ends with 'ba'
        if s.endswith('ba'):
            s = s[:-2] + 'ab'
            steps.append(f"{s} (Rule 2)")
            changed = True
            continue

        # Rule 3: ends with 'ac'
        if s.endswith('ac'):
            s = s[:-2] + 'cb'
            steps.append(f"{s} (Rule 3)")
            changed = True
            continue

        # Rule 4: starts with 'bb'
        if s.startswith('bb'):
            s = s[0] + s[2:]
            steps.append(f"{s} (Rule 4)")
            changed = True
            continue

        # Rule 6: starts with 'ca'
        if s.startswith('ca'):
            s = 'bb' + s[2:] + 'c'
            steps.append(f"{s} (Rule 6)")
            changed = True
            continue

        if not changed:
            break

    print("\nFinal Steps:")
    for i, step in enumerate(steps):
        print(f"Step {i}: {step}")
    print("\nFinal string:", s)

# Test with the string
s = "baccbabbcc"
apply_operations(s)