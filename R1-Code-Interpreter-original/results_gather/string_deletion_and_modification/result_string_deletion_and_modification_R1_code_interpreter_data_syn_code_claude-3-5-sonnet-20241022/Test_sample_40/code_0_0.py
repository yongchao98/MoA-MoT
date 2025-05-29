def apply_rules(s):
    steps = []
    steps.append(f"Initial: {s}")
    
    while True:
        original = s
        changed = False
        
        # Rule 1: Replace first 'abc' with 'cab'
        if 'abc' in s:
            idx = s.index('abc')
            s = s[:idx] + 'cab' + s[idx+3:]
            steps.append(f"Rule 1: {s}")
            changed = True
            continue
            
        # Rule 2: Remove first 'ca' after first character
        if 'ca' in s[1:]:
            idx = s[1:].index('ca') + 1
            s = s[:idx] + s[idx+2:]
            steps.append(f"Rule 2: {s}")
            changed = True
            continue
            
        # Rule 3: Remove last char if ends with 'ca'
        if s.endswith('ca'):
            s = s[:-1]
            steps.append(f"Rule 3: {s}")
            changed = True
            continue
            
        # Rule 4: Remove second char if starts with 'bb'
        if s.startswith('bb'):
            s = s[0] + s[2:]
            steps.append(f"Rule 4: {s}")
            changed = True
            continue
            
        # Rule 5: Replace prefix 'ca' with 'bb' and append 'c'
        if s.startswith('ca'):
            s = 'bb' + s[2:] + 'c'
            steps.append(f"Rule 5: {s}")
            changed = True
            continue
            
        # Rule 6: Replace suffix 'ac' with 'cb'
        if s.endswith('ac'):
            s = s[:-2] + 'cb'
            steps.append(f"Rule 6: {s}")
            changed = True
            continue
            
        if not changed:
            break
    
    for step in steps:
        print(step)
    print("\nFinal string:", s)

# Test with the given string
apply_rules("bbaacaabcbbba")