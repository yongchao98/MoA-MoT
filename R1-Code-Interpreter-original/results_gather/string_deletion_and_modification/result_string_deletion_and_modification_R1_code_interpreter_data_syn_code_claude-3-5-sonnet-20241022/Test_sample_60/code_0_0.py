def apply_rules(s):
    # Returns (new_string, whether_changed)
    
    # Rule 1: Delete first 'bca'
    if 'bca' in s:
        pos = s.find('bca')
        return (s[:pos] + s[pos+3:], True)
    
    # Rule 2: Remove 'ca' after first character
    if len(s) > 1 and 'ca' in s[1:]:
        pos = s[1:].find('ca') + 1  # Add 1 because we started from index 1
        return (s[:pos] + s[pos+2:], True)
    
    # Rule 3: Replace 'aa' at end with 'cc'
    if s.endswith('aa'):
        return (s[:-2] + 'cc', True)
    
    # Rule 4: Replace 'ab' prefix with 'ca'
    if s.startswith('ab'):
        return ('ca' + s[2:], True)
    
    # Rule 5: Replace 'ac' suffix with 'cb'
    if s.endswith('ac'):
        return (s[:-2] + 'cb', True)
    
    # Rule 6: Remove 'cc' prefix
    if s.startswith('cc'):
        return (s[2:], True)
    
    # No rules applied
    return (s, False)

def solve(initial_string):
    current = initial_string
    steps = [current]
    
    while True:
        new_string, changed = apply_rules(current)
        if not changed:
            break
        current = new_string
        steps.append(current)
    
    print("Steps:")
    for i, s in enumerate(steps):
        print(f"Step {i}: {s}")
    return current

# Solve for the given string
initial = "acabacaccbccaab"
result = solve(initial)
print(f"\nFinal result: {result}")