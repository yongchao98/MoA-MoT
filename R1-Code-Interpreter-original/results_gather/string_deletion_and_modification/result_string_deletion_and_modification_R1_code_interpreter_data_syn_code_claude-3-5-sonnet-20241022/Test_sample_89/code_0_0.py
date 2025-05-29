def apply_rules(s):
    steps = []
    steps.append(f"Start: {s}")
    
    while True:
        original = s
        changed = False
        
        # Rule 1: 'ca' not at start
        if s[1:].find('ca') != -1:
            pos = s[1:].find('ca') + 1
            s = s[:pos] + s[pos+2:]
            steps.append(f"Rule 1: {s}")
            changed = True
            continue
            
        # Rule 2: suffix 'bb'
        if s.endswith('bb'):
            s = s[:-2]
            steps.append(f"Rule 2: {s}")
            changed = True
            continue
            
        # Rule 3: ends with 'aa'
        if s.endswith('aa'):
            s = s[:-2] + 'cc'
            steps.append(f"Rule 3: {s}")
            changed = True
            continue
            
        # Rule 4: starts with 'ac'
        if s.startswith('ac'):
            s = 'zz' + s[2:]
            steps.append(f"Rule 4: {s}")
            changed = True
            continue
            
        # Rule 5: contains 'bca'
        if 'bca' in s:
            pos = s.find('bca')
            s = s[:pos] + s[pos+3:]
            steps.append(f"Rule 5: {s}")
            changed = True
            continue
            
        # Rule 6: starts with 'cc'
        if s.startswith('cc'):
            s = s[2:]
            steps.append(f"Rule 6: {s}")
            changed = True
            continue
            
        if not changed:
            break
    
    for step in steps:
        print(step)
    print("\nFinal string:", s)

# Test with the given string
apply_rules("acbbabbcabccacaccc")