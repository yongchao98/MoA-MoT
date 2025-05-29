def apply_operations(s):
    def op1(s):  # Remove middle character if length > 15
        if len(s) > 15:
            mid = len(s) // 2
            return s[:mid] + s[mid+1:], True
        return s, False

    def op2(s):  # Remove 'ca' not at start
        if 'ca' in s[1:]:
            idx = s[1:].index('ca') + 1
            return s[:idx] + s[idx+2:], True
        return s, False

    def op3(s):  # Replace ending 'ba' with 'ab'
        if s.endswith('ba'):
            return s[:-2] + 'ab', True
        return s, False

    def op4(s):  # Replace 'abc' with 'cab'
        if 'abc' in s:
            return s.replace('abc', 'cab', 1), True
        return s, False

    def op5(s):  # Replace 'acb' with 'bca'
        if 'acb' in s:
            return s.replace('acb', 'bca', 1), True
        return s, False

    def op6(s):  # Remove starting 'cc'
        if s.startswith('cc'):
            return s[2:], True
        return s, False

    current = s
    while True:
        original = current
        
        # Apply operations in order
        for op in [op1, op2, op3, op4, op5, op6]:
            current, changed = op(current)
            if changed:
                print(f"Applied operation, new string: {current}")
                break
        
        # If no operation was applied, we're done
        if current == original:
            break
    
    return current

initial = "cacbaaabcac"
result = apply_operations(initial)
print(f"\nFinal result: {result}")