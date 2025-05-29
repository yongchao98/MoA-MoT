def apply_operations(s):
    def op1(s):  # 'ca' not at start
        if len(s) > 2:
            pos = s[1:].find('ca')
            if pos != -1:
                return s[:pos+1] + s[pos+3:], True
        return s, False

    def op2(s):  # ends with 'ba'
        if s.endswith('ba'):
            return s[:-2] + 'ab', True
        return s, False

    def op3(s):  # ends with 'ac'
        if s.endswith('ac'):
            return s[:-2] + 'cb', True
        return s, False

    def op4(s):  # starts with 'bb'
        if s.startswith('bb'):
            return s[0] + s[2:], True
        return s, False

    def op5(s):  # contains 'bca'
        pos = s.find('bca')
        if pos != -1:
            return s[:pos] + s[pos+3:], True
        return s, False

    def op6(s):  # starts with 'ca'
        if s.startswith('ca'):
            return 'bb' + s[2:] + 'c', True
        return s, False

    operations = [op1, op2, op3, op4, op5, op6]
    
    steps = []
    steps.append(s)
    
    while True:
        any_change = False
        for op in operations:
            new_s, changed = op(s)
            if changed:
                s = new_s
                steps.append(s)
                any_change = True
                break
        if not any_change:
            break
    
    print("Initial string:", steps[0])
    print("\nSteps:")
    for i, step in enumerate(steps):
        print(f"Step {i}: {step}")
    print("\nFinal string:", s)

# Test with the given string
s = "cabacaccbacabbcacca"
apply_operations(s)