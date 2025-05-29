def apply_rules(s):
    # Returns (new_string, changed)
    
    # Rule 1: 'ca' not at start
    if 'ca' in s[1:]:
        pos = s[1:].index('ca') + 1
        return s[:pos] + s[pos+2:], True
    
    # Rule 2: prefix 'bc'
    if s.startswith('bc'):
        return s[2:] + 'aa', True
    
    # Rule 3: starts with 'bb'
    if s.startswith('bb'):
        return s[0] + s[2:], True
    
    # Rule 4: prefix 'ca'
    if s.startswith('ca'):
        return 'bb' + s[2:] + 'c', True
    
    # Rule 5: suffix 'cc'
    if s.endswith('cc'):
        return 'a' + s[:-2] + 'b', True
    
    # Rule 6: prefix 'cb'
    if s.startswith('cb'):
        return 'aa' + s[2:-1], True
    
    return s, False

# Initial string
s = "bbcbccacaacbbab"
steps = [s]

while True:
    new_s, changed = apply_rules(s)
    if not changed:
        break
    s = new_s
    steps.append(s)

# Print each step
for i, step in enumerate(steps):
    print(f"Step {i}: {step}")
print("\nFinal string:", s)