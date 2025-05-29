def apply_rules(s):
    # Returns (new_string, changed)
    
    # Rule 1: If prefix is 'ca', replace with 'bb' and append 'c'
    if s.startswith('ca'):
        return s.replace('ca', 'bb', 1) + 'c', True
    
    # Rule 2: If contains 'ca' (not at start), remove first occurrence after first char
    if 'ca' in s[1:]:
        pos = s[1:].find('ca') + 1
        return s[:pos] + s[pos+2:], True
    
    # Rule 3: If length > 15, remove middle character
    if len(s) > 15:
        mid = len(s) // 2
        return s[:mid] + s[mid+1:], True
    
    # Rule 4: If ends with 'aa', replace with 'cc'
    if s.endswith('aa'):
        return s[:-2] + 'cc', True
    
    # Rule 5: If prefix is 'bc', delete first two chars and append 'aa'
    if s.startswith('bc'):
        return s[2:] + 'aa', True
    
    # Rule 6: If starts with 'bb', remove second character
    if s.startswith('bb'):
        return s[0] + s[2:], True
    
    return s, False

def process_string(s):
    steps = [s]
    while True:
        s, changed = apply_rules(s)
        if changed:
            steps.append(s)
        else:
            break
    
    print("Steps:")
    for i, step in enumerate(steps):
        print(f"Step {i}: {step}")
    return s

initial = "bcaccabbbcb"
final = process_string(initial)
print(f"\nFinal result: {final}")