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

    current = s
    steps = []
    steps.append(f"Step 0: {current}")
    step_count = 1
    
    while True:
        original = current
        # Try each operation in order
        for op_num, op in enumerate([op1, op2, op3, op4, op5, op6], 1):
            new_s, changed = op(current)
            if changed:
                current = new_s
                steps.append(f"Step {step_count}: {current} (Rule {op_num})")
                step_count += 1
                break
        # If no operation changed the string, we're done
        if original == current:
            break
    
    for step in steps:
        print(step)
    return current

initial = "cbbcabaccabbcbab"
result = apply_operations(initial)
print(f"\nFinal result: {result}")