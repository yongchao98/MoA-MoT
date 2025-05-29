def continue_operations(s):
    def operation1(s):  # ends with 'ba' -> 'ab'
        if s.endswith('ba'):
            return s[:-2] + 'ab', True
        return s, False

    def operation2(s):  # 'acb' -> 'bca'
        if 'acb' in s:
            idx = s.index('acb')
            return s[:idx] + 'bca' + s[idx+3:], True
        return s, False

    def operation3(s):  # delete 'bca'
        if 'bca' in s:
            idx = s.index('bca')
            return s[:idx] + s[idx+3:], True
        return s, False

    def operation4(s):  # starts with 'ab' -> 'ca'
        if s.startswith('ab'):
            return 'ca' + s[2:], True
        return s, False

    def operation5(s):  # starts with 'bb' -> remove second char
        if s.startswith('bb'):
            return s[0] + s[2:], True
        return s, False

    def operation6(s):  # starts with 'ac' -> 'zz'
        if s.startswith('ac'):
            return 'zz' + s[2:], True
        return s, False

    current = s
    steps = []
    steps.append(current)
    
    while True:
        original = current
        
        # Apply operations in order
        for op in [operation1, operation2, operation3, operation4, operation5, operation6]:
            current, changed = op(current)
            if changed:
                steps.append(current)
                break
        
        # If no operation changed the string, we're done
        if current == original:
            break

    print("Starting string:", s)
    print("\nSteps:")
    for i, step in enumerate(steps):
        print(f"Step {i}: {step}")
    print("\nFinal string:", current)

# Continue from where we left off
s = "caacababbc"
continue_operations(s)