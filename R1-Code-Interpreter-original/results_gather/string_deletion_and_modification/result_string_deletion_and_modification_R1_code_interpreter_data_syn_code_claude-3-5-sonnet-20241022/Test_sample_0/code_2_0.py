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
    
    iteration = 0
    max_iterations = 100  # Safety limit
    
    while iteration < max_iterations:
        old = current
        for i, op in enumerate(operations, 1):
            new_s = op(current)
            if new_s != current:
                print(f"Iteration {iteration}: Rule {i} applies to '{current}' -> '{new_s}'")
                current = new_s
                steps.append(current)
                break
        if old == current:  # No changes were made in this iteration
            break
        iteration += 1
    
    return current

initial = "bbccbbbbacbabbba"
result = apply_operations(initial)
print(f"\nFinal result: {result}")