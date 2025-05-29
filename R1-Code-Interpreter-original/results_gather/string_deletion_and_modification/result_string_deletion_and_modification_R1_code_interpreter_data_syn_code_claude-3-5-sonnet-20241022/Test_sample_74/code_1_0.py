def apply_rules(s):
    # Returns (new_string, changed, rule_applied)
    
    # Rule 1: Replace first 'abc' with 'cab'
    if 'abc' in s:
        pos = s.find('abc')
        return (s[:pos] + 'cab' + s[pos+3:], True, 1)
    
    # Rule 2: Delete first 'bca'
    if 'bca' in s:
        pos = s.find('bca')
        return (s[:pos] + s[pos+3:], True, 2)
    
    # Rule 3: If even number of b's, append 'ab'
    b_count = s.count('b')
    if b_count > 0 and b_count % 2 == 0:
        return (s + 'ab', True, 3)
    
    # Rule 4: If starts with 'cc', remove first two chars
    if s.startswith('cc'):
        return (s[2:], True, 4)
    
    # Rule 5: If starts with 'ca', replace with 'bb' and append 'c'
    if s.startswith('ca'):
        return ('bb' + s[2:] + 'c', True, 5)
    
    # Rule 6: If starts with 'bc', delete first two chars and append 'aa'
    if s.startswith('bc'):
        return (s[2:] + 'aa', True, 6)
    
    return (s, False, 0)

# Initial string
s = "ccaabbcaccbbcbbbaaa"
steps = []
steps.append(("Initial", s))

while True:
    new_s, changed, rule = apply_rules(s)
    if not changed:
        break
    s = new_s
    steps.append((f"Rule {rule}", s))

# Print all steps and final result
for i, (rule, step) in enumerate(steps):
    print(f"Step {i}: {step} ({rule})")
print("\nFinal string:", s)