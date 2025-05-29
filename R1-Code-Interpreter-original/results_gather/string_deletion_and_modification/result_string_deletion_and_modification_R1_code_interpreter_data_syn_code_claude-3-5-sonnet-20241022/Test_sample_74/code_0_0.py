def apply_rules(s):
    # Returns (new_string, changed)
    
    # Rule 1: Replace first 'abc' with 'cab'
    if 'abc' in s:
        pos = s.find('abc')
        return (s[:pos] + 'cab' + s[pos+3:], True)
    
    # Rule 2: Delete first 'bca'
    if 'bca' in s:
        pos = s.find('bca')
        return (s[:pos] + s[pos+3:], True)
    
    # Rule 3: If even number of b's, append 'ab'
    b_count = s.count('b')
    if b_count > 0 and b_count % 2 == 0:
        return (s + 'ab', True)
    
    # Rule 4: If starts with 'cc', remove first two chars
    if s.startswith('cc'):
        return (s[2:], True)
    
    # Rule 5: If starts with 'ca', replace with 'bb' and append 'c'
    if s.startswith('ca'):
        return ('bb' + s[2:] + 'c', True)
    
    # Rule 6: If starts with 'bc', delete first two chars and append 'aa'
    if s.startswith('bc'):
        return (s[2:] + 'aa', True)
    
    return (s, False)

# Initial string
s = "ccaabbcaccbbcbbbaaa"
steps = []
steps.append(s)

while True:
    new_s, changed = apply_rules(s)
    if not changed:
        break
    s = new_s
    steps.append(s)

# Print all steps and final result
for i, step in enumerate(steps):
    print(f"Step {i}: {step}")
print("\nFinal string:", s)