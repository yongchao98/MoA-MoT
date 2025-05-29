def apply_rules(s):
    # Returns (new_string, changed)
    # Rule 1: If prefix is 'aa', remove first character
    if s.startswith('aa'):
        return s[1:], True
    
    # Rule 2: If ends with 'aa', replace with 'cc'
    if s.endswith('aa'):
        return s[:-2] + 'cc', True
    
    # Rule 3: If prefix is 'ca', replace with 'bb' and append 'c'
    if s.startswith('ca'):
        return 'bb' + s[2:] + 'c', True
    
    # Rule 4: If suffix is 'cc', replace with 'b' and prepend 'a'
    if s.endswith('cc'):
        return 'a' + s[:-2] + 'b', True
    
    # Rule 5: If suffix is 'bb', delete last two characters
    if s.endswith('bb'):
        return s[:-2], True
    
    # Rule 6: If ends with 'ca', remove last character
    if s.endswith('ca'):
        return s[:-1], True
    
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

# Test with the given string
initial = "aabccbccbbccbb"
final = process_string(initial)
print(f"\nFinal result: \"{final}\"")