def apply_rules(s):
    def count_b(s):
        return s.count('b')
    
    changes = True
    steps = []
    steps.append(f"Start: {s}")
    
    while changes:
        changes = False
        
        # Rule 1: If starts with 'bb', remove second character
        if s.startswith('bb'):
            s = s[0] + s[2:]
            steps.append(f"Rule 1: {s}")
            changes = True
            continue
            
        # Rule 2: If prefix is 'ab', replace with 'ca'
        if s.startswith('ab'):
            s = 'ca' + s[2:]
            steps.append(f"Rule 2: {s}")
            changes = True
            continue
            
        # Rule 3: If contains 'bca', delete first occurrence
        if 'bca' in s:
            idx = s.index('bca')
            s = s[:idx] + s[idx+3:]
            steps.append(f"Rule 3: {s}")
            changes = True
            continue
            
        # Rule 4: If even number of b's (and at least one), append 'ab'
        b_count = count_b(s)
        if b_count > 0 and b_count % 2 == 0:
            s = s + 'ab'
            steps.append(f"Rule 4: {s}")
            changes = True
            continue
            
        # Rule 5: If prefix is 'cb', replace with 'aa' and delete last char
        if s.startswith('cb'):
            s = 'aa' + s[2:]
            s = s[:-1]
            steps.append(f"Rule 5: {s}")
            changes = True
            continue
            
        # Rule 6: If suffix is 'ac', replace with 'cb'
        if s.endswith('ac'):
            s = s[:-2] + 'cb'
            steps.append(f"Rule 6: {s}")
            changes = True
            continue
    
    print("\n".join(steps))
    print("\nFinal string:", s)

# Test with the given string
apply_rules("bbbcaaaabbbaccaab")