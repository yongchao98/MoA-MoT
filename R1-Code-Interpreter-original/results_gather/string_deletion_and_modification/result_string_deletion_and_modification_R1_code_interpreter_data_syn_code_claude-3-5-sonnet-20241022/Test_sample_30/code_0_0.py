def apply_operations(s):
    def op1(s):  # If prefix is 'aa', remove first character
        return s[1:] if s.startswith('aa') else s
    
    def op2(s):  # If contains 'abc', replace first occurrence with 'cab'
        return s.replace('abc', 'cab', 1) if 'abc' in s else s
    
    def op3(s):  # If starts with 'bb', remove second character
        return s[0] + s[2:] if s.startswith('bb') else s
    
    def op4(s):  # If starts with 'cc', remove first two characters
        return s[2:] if s.startswith('cc') else s
    
    def op5(s):  # If ends with 'ba', replace with 'ab'
        return s[:-2] + 'ab' if s.endswith('ba') else s
    
    def op6(s):  # If suffix is 'cc', replace with 'b' and prepend 'a'
        return 'a' + s[:-2] + 'b' if s.endswith('cc') else s

    operations = [op1, op2, op3, op4, op5, op6]
    current = s
    steps = []
    steps.append(current)
    
    while True:
        changed = False
        for op in operations:
            new_s = op(current)
            if new_s != current:
                current = new_s
                steps.append(current)
                changed = True
                break
        if not changed:
            break
    
    print("Initial string:", s)
    print("\nSteps:")
    for i, step in enumerate(steps):
        print(f"Step {i}: {step}")
    print("\nFinal string:", current)

# Test with the given string
apply_operations("aaacbcaccacaacbcc")