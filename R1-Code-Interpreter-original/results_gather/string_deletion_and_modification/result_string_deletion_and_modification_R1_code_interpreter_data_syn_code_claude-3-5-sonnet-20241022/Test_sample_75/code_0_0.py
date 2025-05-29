def apply_operations(s):
    def operation1(s):  # suffix 'bb' -> delete last two
        if s.endswith('bb'):
            return s[:-2], True
        return s, False
    
    def operation2(s):  # 'acb' -> 'bca'
        if 'acb' in s:
            idx = s.index('acb')
            return s[:idx] + 'bca' + s[idx+3:], True
        return s, False
    
    def operation3(s):  # suffix 'aa' -> 'cc'
        if s.endswith('aa'):
            return s[:-2] + 'cc', True
        return s, False
    
    def operation4(s):  # 'bca' -> delete
        if 'bca' in s:
            idx = s.index('bca')
            return s[:idx] + s[idx+3:], True
        return s, False
    
    def operation5(s):  # suffix 'cc' -> prepend 'a' and replace with 'b'
        if s.endswith('cc'):
            return 'a' + s[:-2] + 'b', True
        return s, False
    
    def operation6(s):  # suffix 'ba' -> 'ab'
        if s.endswith('ba'):
            return s[:-2] + 'ab', True
        return s, False

    operations = [operation1, operation2, operation3, operation4, operation5, operation6]
    
    current = s
    steps = []
    steps.append(current)
    
    while True:
        any_change = False
        for op in operations:
            new_s, changed = op(current)
            if changed:
                current = new_s
                steps.append(current)
                any_change = True
                break
        if not any_change:
            break
    
    print("Steps:")
    for i, step in enumerate(steps):
        print(f"Step {i}: {step}")
    
    return current

initial = "acaccacbcbcbabbbacb"
result = apply_operations(initial)
print(f"\nFinal result: {result}")