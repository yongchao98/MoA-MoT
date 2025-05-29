def apply_operations(s):
    def op1(s):  # suffix 'ac' -> 'cb'
        return s[:-2] + 'cb' if s.endswith('ac') else s
    
    def op2(s):  # starts with 'cc' -> remove first two chars
        return s[2:] if s.startswith('cc') else s
    
    def op3(s):  # 'abc' -> 'cab'
        return s.replace('abc', 'cab', 1) if 'abc' in s else s
    
    def op4(s):  # prefix 'bc' -> delete first two and append 'aa'
        return s[2:] + 'aa' if s.startswith('bc') else s
    
    def op5(s):  # 'ca' not at start -> remove first occurrence after first char
        if 'ca' in s[1:]:
            idx = s.find('ca', 1)
            return s[:idx] + s[idx+2:]
        return s
    
    def op6(s):  # suffix 'bb' -> delete last two chars
        return s[:-2] if s.endswith('bb') else s
    
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
    
    print("Steps:")
    for i, step in enumerate(steps):
        print(f"Step {i}: {step}")
    return current

initial = "bbccacbbbbaabccabbba"
result = apply_operations(initial)
print(f"\nFinal result: {result}")