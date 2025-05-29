def apply_rules(s):
    def rule1(s):  # 'ac' suffix -> 'cb'
        if s.endswith('ac'):
            return s[:-2] + 'cb'
        return s

    def rule2(s):  # 'abc' -> 'cab'
        if 'abc' in s:
            i = s.index('abc')
            return s[:i] + 'cab' + s[i+3:]
        return s

    def rule3(s):  # 'cc' suffix -> prepend 'a' and replace with 'b'
        if s.endswith('cc'):
            return 'a' + s[:-2] + 'b'
        return s

    def rule4(s):  # 'ba' suffix -> 'ab'
        if s.endswith('ba'):
            return s[:-2] + 'ab'
        return s

    def rule5(s):  # 'bc' prefix -> delete and append 'aa'
        if s.startswith('bc'):
            return s[2:] + 'aa'
        return s

    def rule6(s):  # 'ca' not at start -> remove first occurrence after start
        if len(s) > 1:
            i = s[1:].find('ca')
            if i != -1:
                i += 1  # adjust index for the slice we took
                return s[:i] + s[i+2:]
        return s

    original = s
    steps = []
    steps.append(f"Start: {s}")
    
    while True:
        changed = False
        
        # Try each rule in order
        new_s = rule1(s)
        if new_s != s:
            steps.append(f"Rule 1: {s} -> {new_s}")
            s = new_s
            changed = True
            continue
            
        new_s = rule2(s)
        if new_s != s:
            steps.append(f"Rule 2: {s} -> {new_s}")
            s = new_s
            changed = True
            continue
            
        new_s = rule3(s)
        if new_s != s:
            steps.append(f"Rule 3: {s} -> {new_s}")
            s = new_s
            changed = True
            continue
            
        new_s = rule4(s)
        if new_s != s:
            steps.append(f"Rule 4: {s} -> {new_s}")
            s = new_s
            changed = True
            continue
            
        new_s = rule5(s)
        if new_s != s:
            steps.append(f"Rule 5: {s} -> {new_s}")
            s = new_s
            changed = True
            continue
            
        new_s = rule6(s)
        if new_s != s:
            steps.append(f"Rule 6: {s} -> {new_s}")
            s = new_s
            changed = True
            continue
        
        if not changed:
            break
    
    for step in steps:
        print(step)
    print(f"\nFinal string: {s}")

# Test with the given string
apply_rules("ccacacaacc")