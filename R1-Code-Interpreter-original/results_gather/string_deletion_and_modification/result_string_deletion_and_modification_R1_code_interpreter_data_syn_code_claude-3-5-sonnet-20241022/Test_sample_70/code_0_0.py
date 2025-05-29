def apply_rules(s):
    # Returns (new_string, changed)
    
    # Rule 1: Remove 'ca' not at start
    if 'ca' in s[1:]:
        pos = s[1:].index('ca') + 1  # Add 1 because we started search from index 1
        return s[:pos] + s[pos+2:], True
    
    # Rule 2: Replace 'acb' with 'bca'
    if 'acb' in s:
        pos = s.index('acb')
        return s[:pos] + 'bca' + s[pos+3:], True
    
    # Rule 3: Replace 'abc' with 'cab'
    if 'abc' in s:
        pos = s.index('abc')
        return s[:pos] + 'cab' + s[pos+3:], True
    
    # Rule 4: If starts with 'bc', delete first two and append 'aa'
    if s.startswith('bc'):
        return s[2:] + 'aa', True
    
    # Rule 5: If starts with 'aa', remove first character
    if s.startswith('aa'):
        return s[1:], True
    
    # Rule 6: If starts with 'cc', remove first two characters
    if s.startswith('cc'):
        return s[2:], True
    
    return s, False

# Initial string
s = "acaabacaccbcbccbcaa"

# Keep track of steps
steps = [s]
changed = True

while changed:
    s, changed = apply_rules(s)
    if changed:
        steps.append(s)

# Print all steps and final result
print("Steps:")
for i, step in enumerate(steps):
    print(f"Step {i}: {step}")
print("\nFinal string:", s)