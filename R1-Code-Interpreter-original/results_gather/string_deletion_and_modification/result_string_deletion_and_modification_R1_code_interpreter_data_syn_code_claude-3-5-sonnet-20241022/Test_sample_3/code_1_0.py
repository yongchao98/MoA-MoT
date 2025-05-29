def apply_rule1(s):
    # If suffix is 'ac', replace with 'cb'
    if s.endswith('ac'):
        return s[:-2] + 'cb'
    return s

def apply_rule2(s):
    # If length > 15, remove middle character
    if len(s) > 15:
        mid = len(s) // 2
        return s[:mid] + s[mid+1:]
    return s

def apply_rule3(s):
    # If ends with 'ca', remove last character
    if s.endswith('ca'):
        return s[:-1]
    return s

def apply_rule4(s):
    # If even number of 'b's (and at least one), append 'ab'
    b_count = s.count('b')
    if b_count > 0 and b_count % 2 == 0:
        return s + 'ab'
    return s

def apply_rule5(s):
    # If contains 'bca', delete first occurrence
    if 'bca' in s:
        idx = s.find('bca')
        return s[:idx] + s[idx+3:]
    return s

def apply_rule6(s):
    # If starts with 'ac', replace first two chars with 'zz'
    if s.startswith('ac'):
        return 'zz' + s[2:]
    return s

def solve(s):
    prev = ""
    current = s
    steps = 1
    
    while prev != current and len(current) > 0:
        prev = current
        print(f"\nStep {steps}:")
        print("Starting:", current)
        
        # Apply each rule in sequence
        for i in range(1, 7):
            old = current
            if i == 1:
                current = apply_rule1(current)
            elif i == 2:
                current = apply_rule2(current)
            elif i == 3:
                current = apply_rule3(current)
            elif i == 4:
                current = apply_rule4(current)
            elif i == 5:
                current = apply_rule5(current)
            else:
                current = apply_rule6(current)
            print(f"After rule {i}: {current}")
            
        steps += 1
    
    return current

# Initial string
initial = "cbaccccccbcabcbcabcb"
result = solve(initial)
print("\nFinal result:", result)