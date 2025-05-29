def apply_rules(s):
    steps = [s]
    while True:
        original = s
        
        # Rule 1: If starts with 'bb', remove second character
        if s.startswith('bb'):
            s = s[0] + s[2:]
            steps.append(s)
            continue
            
        # Rule 2: If length > 15, remove middle character
        if len(s) > 15:
            mid = len(s) // 2
            s = s[:mid] + s[mid+1:]
            steps.append(s)
            continue
            
        # Rule 3: If ends with 'aa', replace with 'cc'
        if s.endswith('aa'):
            s = s[:-2] + 'cc'
            steps.append(s)
            continue
            
        # Rule 4: If starts with 'bc', delete first two and append 'aa'
        if s.startswith('bc'):
            s = s[2:] + 'aa'
            steps.append(s)
            continue
            
        # Rule 5: If even number of 'b's (and at least one), append 'ab'
        b_count = s.count('b')
        if b_count > 0 and b_count % 2 == 0:
            s = s + 'ab'
            steps.append(s)
            continue
            
        # Rule 6: If ends with 'bb', delete last two characters
        if s.endswith('bb'):
            s = s[:-2]
            steps.append(s)
            continue
            
        if s == original:
            break
    
    print("Steps:")
    for i, step in enumerate(steps):
        print(f"Step {i}: {step}")
    print("\nFinal string:", s)

initial = "ababcaccbbbccbaabcbc"
apply_rules(initial)