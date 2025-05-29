def apply_operations(s):
    def op1(s):  # If starts with 'bb', remove second character
        if s.startswith('bb'):
            return s[0] + s[2:], True
        return s, False
    
    def op2(s):  # If contains 'acb', replace first occurrence with 'bca'
        if 'acb' in s:
            idx = s.index('acb')
            return s[:idx] + 'bca' + s[idx+3:], True
        return s, False
    
    def op3(s):  # If suffix is 'bb', delete last two characters
        if s.endswith('bb'):
            return s[:-2], True
        return s, False
    
    def op4(s):  # If suffix is 'cc', replace with 'b' and prepend 'a'
        if s.endswith('cc'):
            return 'a' + s[:-2] + 'b', True
        return s, False
    
    def op5(s):  # If contains 'abc', replace first occurrence with 'cab'
        if 'abc' in s:
            idx = s.index('abc')
            return s[:idx] + 'cab' + s[idx+3:], True
        return s, False
    
    def op6(s):  # If ends with 'ba', replace with 'ab'
        if s.endswith('ba'):
            return s[:-2] + 'ab', True
        return s, False

    operations = [op1, op2, op3, op4, op5, op6]
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

initial = "cbbaacbabccbcbba"
result = apply_operations(initial)
print(f"\nFinal result: {result}")