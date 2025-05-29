def apply_operations(s):
    def op1(s):  # If starts with 'bb', remove second character
        if s.startswith('bb'):
            return s[0] + s[2:], True, "Rule 1: bb* -> b*"
        return s, False, ""
    
    def op2(s):  # If contains 'acb', replace first occurrence with 'bca'
        if 'acb' in s:
            idx = s.index('acb')
            return s[:idx] + 'bca' + s[idx+3:], True, "Rule 2: acb -> bca"
        return s, False, ""
    
    def op3(s):  # If suffix is 'bb', delete last two characters
        if s.endswith('bb'):
            return s[:-2], True, "Rule 3: *bb -> *"
        return s, False, ""
    
    def op4(s):  # If suffix is 'cc', replace with 'b' and prepend 'a'
        if s.endswith('cc'):
            return 'a' + s[:-2] + 'b', True, "Rule 4: *cc -> a*b"
        return s, False, ""
    
    def op5(s):  # If contains 'abc', replace first occurrence with 'cab'
        if 'abc' in s:
            idx = s.index('abc')
            return s[:idx] + 'cab' + s[idx+3:], True, "Rule 5: abc -> cab"
        return s, False, ""
    
    def op6(s):  # If ends with 'ba', replace with 'ab'
        if s.endswith('ba'):
            return s[:-2] + 'ab', True, "Rule 6: *ba -> *ab"
        return s, False, ""

    operations = [op1, op2, op3, op4, op5, op6]
    current = "cbbcabaccabbcbab"
    steps = []
    steps.append(f"Start: {current}")
    
    while True:
        any_change = False
        for op in operations:
            new_s, changed, rule = op(current)
            if changed:
                current = new_s
                steps.append(f"{current} ({rule})")
                any_change = True
                break
        if not any_change:
            break
    
    for step in steps:
        print(step)
    return current

result = apply_operations("cbbcabaccabbcbab")
print(f"\nFinal result: {result}")