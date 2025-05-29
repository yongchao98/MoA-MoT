def operation1(s):
    # If contains 'ca' (not at start), remove first occurrence after first char
    pos = s[1:].find('ca')
    if pos != -1:
        return s[:pos+1] + s[pos+3:]
    return s

def operation2(s):
    # If even number of b's (at least one), append 'ab'
    b_count = s.count('b')
    if b_count > 0 and b_count % 2 == 0:
        return s + 'ab'
    return s

def operation3(s):
    # If contains 'bca', delete first occurrence
    pos = s.find('bca')
    if pos != -1:
        return s[:pos] + s[pos+3:]
    return s

def operation4(s):
    # If prefix is 'cb', replace with 'aa' and delete last char
    if s.startswith('cb'):
        return 'aa' + s[2:-1]
    return s

def operation5(s):
    # If prefix is 'ab', replace with 'ca'
    if s.startswith('ab'):
        return 'ca' + s[2:]
    return s

def operation6(s):
    # If suffix is 'cc', replace with 'b' and prepend 'a'
    if s.endswith('cc'):
        return 'a' + s[:-2] + 'b'
    return s

def apply_operations(s):
    operations = [operation1, operation2, operation3, operation4, operation5, operation6]
    changed = True
    steps = []
    
    while changed:
        changed = False
        for op in operations:
            new_s = op(s)
            if new_s != s:
                s = new_s
                steps.append(s)
                changed = True
                break
    
    print("Initial string:", "acaccbcbaacccbcccaac")
    print("\nSteps:")
    for i, step in enumerate(steps, 1):
        print(f"Step {i}: {step}")
    print("\nFinal string:", s)

# Run the operations
initial = "acaccbcbaacccbcccaac"
apply_operations(initial)