def apply_rules(s):
    # Returns (new_string, changed)
    
    # Rule 1: If suffix is 'bb', delete last two characters
    if s.endswith('bb'):
        return s[:-2], True
    
    # Rule 2: If suffix is 'ac', replace with 'cb'
    if s.endswith('ac'):
        return s[:-2] + 'cb', True
    
    # Rule 3: If contains 'acb', replace first occurrence with 'bca'
    if 'acb' in s:
        pos = s.find('acb')
        return s[:pos] + 'bca' + s[pos+3:], True
    
    # Rule 4: If prefix is 'bc', delete first two and append 'aa'
    if s.startswith('bc'):
        return s[2:] + 'aa', True
    
    # Rule 5: If prefix is 'ab', replace with 'ca'
    if s.startswith('ab'):
        return 'ca' + s[2:], True
    
    # Rule 6: If prefix is 'ca', replace with 'bb' and append 'c'
    if s.startswith('ca'):
        return 'bb' + s[2:] + 'c', True
    
    return s, False

# Initial string
s = "bccaaaacaa"
steps = []
steps.append(s)

while True:
    new_s, changed = apply_rules(s)
    if not changed:
        break
    s = new_s
    steps.append(s)

# Print all steps and final result
print("Steps:")
for i, step in enumerate(steps):
    print(f"Step {i}: {step}")
print("\nFinal string:", s)