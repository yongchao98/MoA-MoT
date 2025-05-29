def apply_rules(s):
    steps = [s]
    while True:
        changed = False
        
        # Store original for comparison
        original = s
        
        # Print current state and counts for debugging
        print(f"\nCurrent string: {s}")
        print(f"Length: {len(s)}")
        print(f"Number of b's: {s.count('b')}")
        print(f"Starts with bb: {s.startswith('bb')}")
        print(f"Ends with aa: {s.endswith('aa')}")
        print(f"Starts with bc: {s.startswith('bc')}")
        print(f"Ends with bb: {s.endswith('bb')}")
        
        # Rule 1: If starts with 'bb', remove second character
        if s.startswith('bb'):
            s = s[0] + s[2:]
            steps.append(s)
            print("Applied Rule 1")
            continue
            
        # Rule 2: If length > 15, remove middle character
        if len(s) > 15:
            mid = len(s) // 2
            s = s[:mid] + s[mid+1:]
            steps.append(s)
            print("Applied Rule 2")
            continue
            
        # Rule 3: If ends with 'aa', replace with 'cc'
        if s.endswith('aa'):
            s = s[:-2] + 'cc'
            steps.append(s)
            print("Applied Rule 3")
            continue
            
        # Rule 4: If starts with 'bc', delete first two and append 'aa'
        if s.startswith('bc'):
            s = s[2:] + 'aa'
            steps.append(s)
            print("Applied Rule 4")
            continue
            
        # Rule 5: If even number of 'b's (and at least one), append 'ab'
        b_count = s.count('b')
        if b_count > 0 and b_count % 2 == 0:
            s = s + 'ab'
            steps.append(s)
            print("Applied Rule 5")
            continue
            
        # Rule 6: If ends with 'bb', delete last two characters
        if s.endswith('bb'):
            s = s[:-2]
            steps.append(s)
            print("Applied Rule 6")
            continue
            
        if s == original:
            break
    
    print("\nFinal Steps:")
    for i, step in enumerate(steps):
        print(f"Step {i}: {step}")
    print("\nFinal string:", s)

initial = "ababcaccbbbccbaabcbc"
apply_rules(initial)