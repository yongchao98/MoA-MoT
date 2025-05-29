def apply_operations(s):
    def operation1(s):  # 'ca' not at start
        if len(s) > 2:
            pos = s[1:].find('ca')
            if pos != -1:
                return s[:pos+1] + s[pos+3:], True
        return s, False

    def operation2(s):  # starts with 'ac'
        if s.startswith('ac'):
            return 'zz' + s[2:], True
        return s, False

    def operation3(s):  # ends with 'ba'
        if s.endswith('ba'):
            return s[:-2] + 'ab', True
        return s, False

    def operation4(s):  # starts with 'cb'
        if s.startswith('cb'):
            return 'aa' + s[2:-1], True
        return s, False

    def operation5(s):  # starts with 'cc'
        if s.startswith('cc'):
            return s[2:], True
        return s, False

    def operation6(s):  # starts with 'bb'
        if s.startswith('bb'):
            return s[0] + s[2:], True
        return s, False

    current = s
    steps = []
    steps.append(current)
    
    while True:
        changed = False
        
        # Try each operation in order
        s, change = operation1(current)
        if change:
            current = s
            steps.append(f"After op1: {current}")
            continue
            
        s, change = operation2(current)
        if change:
            current = s
            steps.append(f"After op2: {current}")
            continue
            
        s, change = operation3(current)
        if change:
            current = s
            steps.append(f"After op3: {current}")
            continue
            
        s, change = operation4(current)
        if change:
            current = s
            steps.append(f"After op4: {current}")
            continue
            
        s, change = operation5(current)
        if change:
            current = s
            steps.append(f"After op5: {current}")
            continue
            
        s, change = operation6(current)
        if change:
            current = s
            steps.append(f"After op6: {current}")
            continue
            
        if not change:
            break
    
    print("\n".join(steps))
    print("\nFinal string:", current)

# Test with the given string
apply_operations("cbbcbcacccaabbcbcbca")