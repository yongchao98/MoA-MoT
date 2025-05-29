def apply_operations(s):
    def operation1(s):  # 'ca' not at start
        if len(s) > 1:
            pos = s[1:].find('ca')
            if pos != -1:
                return s[:pos+1] + s[pos+3:], True
        return s, False

    def operation2(s):  # prefix 'ca'
        if s.startswith('ca'):
            return 'bb' + s[2:] + 'c', True
        return s, False

    def operation3(s):  # contains 'abc'
        if 'abc' in s:
            pos = s.find('abc')
            return s[:pos] + 'cab' + s[pos+3:], True
        return s, False

    def operation4(s):  # suffix 'ac'
        if s.endswith('ac'):
            return s[:-2] + 'cb', True
        return s, False

    def operation5(s):  # ends with 'ba'
        if s.endswith('ba'):
            return s[:-2] + 'ab', True
        return s, False

    def operation6(s):  # starts with 'cc'
        if s.startswith('cc'):
            return s[2:], True
        return s, False

    current = s
    steps = []
    steps.append(current)
    
    while True:
        changed = False
        
        # Try each operation in order
        for op in [operation1, operation2, operation3, operation4, operation5, operation6]:
            new_s, was_changed = op(current)
            if was_changed:
                current = new_s
                steps.append(current)
                changed = True
                break
        
        if not changed:
            break
    
    print("Steps:")
    for i, step in enumerate(steps):
        print(f"Step {i}: {step}")
    
    return current

initial = "acccbbaaaacacacbbcba"
result = apply_operations(initial)
print("\nFinal result:", result)