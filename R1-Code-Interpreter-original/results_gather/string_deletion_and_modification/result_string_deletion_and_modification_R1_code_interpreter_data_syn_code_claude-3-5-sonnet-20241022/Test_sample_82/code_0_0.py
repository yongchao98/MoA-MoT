def operation1(s):
    # If ends with 'ac', replace with 'cb'
    if s.endswith('ac'):
        return s[:-2] + 'cb', True
    return s, False

def operation2(s):
    # If even number of b's (and at least one), append 'ab'
    b_count = s.count('b')
    if b_count > 0 and b_count % 2 == 0:
        return s + 'ab', True
    return s, False

def operation3(s):
    # If contains 'ca' (not at start), remove first occurrence after first char
    if 'ca' in s[1:]:
        idx = s.find('ca', 1)
        return s[:idx] + s[idx+2:], True
    return s, False

def operation4(s):
    # If length > 15, remove middle character
    if len(s) > 15:
        mid = len(s) // 2
        return s[:mid] + s[mid+1:], True
    return s, False

def operation5(s):
    # If contains 'acb', replace first occurrence with 'bca'
    if 'acb' in s:
        idx = s.find('acb')
        return s[:idx] + 'bca' + s[idx+3:], True
    return s, False

def operation6(s):
    # If contains 'bca', delete first occurrence
    if 'bca' in s:
        idx = s.find('bca')
        return s[:idx] + s[idx+3:], True
    return s, False

def process_string(s):
    operations = [operation1, operation2, operation3, operation4, operation5, operation6]
    steps = []
    steps.append(s)
    
    while True:
        original = s
        for op in operations:
            s, changed = op(s)
            if changed:
                steps.append(s)
                break
        if original == s:  # No changes made in full cycle
            break
    
    print("Initial string:", steps[0])
    for i, step in enumerate(steps[1:], 1):
        print(f"Step {i}: {step}")
    print("\nFinal string:", s)

# Process the given string
process_string("cbbcaaaaaaccbabaabbc")