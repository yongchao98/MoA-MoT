def apply_rules(s):
    original = s
    
    # Rule 1: Replace first 'acb' with 'bca'
    if 'acb' in s:
        pos = s.find('acb')
        s = s[:pos] + 'bca' + s[pos+3:]
        return s
    
    # Rule 2: If starts with 'ab', replace with 'ca'
    if s.startswith('ab'):
        s = 'ca' + s[2:]
        return s
    
    # Rule 3: If ends with 'ba', replace with 'ab'
    if s.endswith('ba'):
        s = s[:-2] + 'ab'
        return s
    
    # Rule 4: If ends with 'ca', remove last character
    if s.endswith('ca'):
        s = s[:-1]
        return s
    
    # Rule 5: If starts with 'ca', replace with 'bb' and append 'c'
    if s.startswith('ca'):
        s = 'bb' + s[2:] + 'c'
        return s
    
    # Rule 6: If contains 'ca' (not at start), remove first such occurrence
    if 'ca' in s[1:]:
        pos = s[1:].find('ca') + 1
        s = s[:pos] + s[pos+2:]
        return s
    
    return s

# Initial string
s = "cacbcbcbac"
steps = [s]

# Keep applying rules until no changes occur
while True:
    new_s = apply_rules(s)
    if new_s == s:
        break
    s = new_s
    steps.append(s)

# Print each step
for i, step in enumerate(steps):
    print(f"Step {i}: {step}")
print("\nFinal string:", s)