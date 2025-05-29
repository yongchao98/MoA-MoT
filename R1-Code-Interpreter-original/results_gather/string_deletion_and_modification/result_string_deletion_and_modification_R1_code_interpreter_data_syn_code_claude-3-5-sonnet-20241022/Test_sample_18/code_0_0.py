def apply_operations(s):
    def operation1(s):  # starts with 'cc'
        return s[2:] if s.startswith('cc') else s
    
    def operation2(s):  # prefix is 'aa'
        return s[1:] if s.startswith('aa') else s
    
    def operation3(s):  # contains 'bca'
        if 'bca' in s:
            i = s.index('bca')
            return s[:i] + s[i+3:]
        return s
    
    def operation4(s):  # suffix is 'bb'
        return s[:-2] if s.endswith('bb') else s
    
    def operation5(s):  # contains 'ca' not at start
        if len(s) > 1 and 'ca' in s[1:]:
            i = s.index('ca', 1)
            return s[:i] + s[i+2:]
        return s
    
    def operation6(s):  # ends with 'aa'
        return s[:-2] + 'cc' if s.endswith('aa') else s

    current = s
    steps = []
    steps.append(f"Start: {current}")
    
    while True:
        original = current
        
        # Apply operations in order
        current = operation1(current)
        if current != original:
            steps.append(f"Op1: {current}")
            continue
            
        current = operation2(current)
        if current != original:
            steps.append(f"Op2: {current}")
            continue
            
        current = operation3(current)
        if current != original:
            steps.append(f"Op3: {current}")
            continue
            
        current = operation4(current)
        if current != original:
            steps.append(f"Op4: {current}")
            continue
            
        current = operation5(current)
        if current != original:
            steps.append(f"Op5: {current}")
            continue
            
        current = operation6(current)
        if current != original:
            steps.append(f"Op6: {current}")
            continue
            
        break
    
    for step in steps:
        print(step)
    print("\nFinal result:", current)

s = "bbbbbacacbccaaaab"
apply_operations(s)