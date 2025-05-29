def apply_rules(s):
    # Returns (new_string, changed)
    
    # Rule 1: If prefix is 'ab', replace with 'ca'
    if s.startswith('ab'):
        return 'ca' + s[2:], True
    
    # Rule 2: If suffix is 'bb', delete last two characters
    if s.endswith('bb'):
        return s[:-2], True
    
    # Rule 3: If ends with 'ba', replace with 'ab'
    if s.endswith('ba'):
        return s[:-2] + 'ab', True
    
    # Rule 4: If length > 15, remove middle character
    if len(s) > 15:
        mid = len(s) // 2
        return s[:mid] + s[mid+1:], True
    
    # Rule 5: If suffix is 'ac', replace with 'cb'
    if s.endswith('ac'):
        return s[:-2] + 'cb', True
    
    # Rule 6: If prefix is 'ca', replace with 'bb' and append 'c'
    if s.startswith('ca'):
        return 'bb' + s[2:] + 'c', True
    
    return s, False

def process_string(s):
    steps = [s]
    while True:
        new_s, changed = apply_rules(s)
        if not changed:
            break
        s = new_s
        steps.append(s)
    
    print("Steps:")
    for i, step in enumerate(steps):
        print(f"Step {i}: {step}")
    return s

initial = "bcbcabcababcbacacbab"
final = process_string(initial)
print(f"\nFinal result: \"{final}\"")