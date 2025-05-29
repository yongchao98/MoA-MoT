def apply_rules(s):
    # Returns (new_string, changed)
    
    # Rule 1: If prefix is 'bc', delete first two and append 'aa'
    if s.startswith('bc'):
        return s[2:] + 'aa', True
    
    # Rule 2: If ends with 'aa', replace with 'cc'
    if s.endswith('aa'):
        return s[:-2] + 'cc', True
    
    # Rule 3: If prefix is 'aa', remove first character
    if s.startswith('aa'):
        return s[1:], True
    
    # Rule 4: If contains 'ca' (not at start), remove first occurrence after first char
    ca_pos = s[1:].find('ca')
    if ca_pos != -1:
        ca_pos += 1  # Adjust position to account for full string
        return s[:ca_pos] + s[ca_pos+2:], True
    
    # Rule 5: If prefix is 'cb', replace with 'aa' and delete last char
    if s.startswith('cb'):
        return 'aa' + s[2:-1], True
    
    # Rule 6: If contains 'acb', replace first occurrence with 'bca'
    acb_pos = s.find('acb')
    if acb_pos != -1:
        return s[:acb_pos] + 'bca' + s[acb_pos+3:], True
    
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

# Process the given string
initial = "cbcbaabacc"
final = process_string(initial)
print(f"\nFinal result: {final}")