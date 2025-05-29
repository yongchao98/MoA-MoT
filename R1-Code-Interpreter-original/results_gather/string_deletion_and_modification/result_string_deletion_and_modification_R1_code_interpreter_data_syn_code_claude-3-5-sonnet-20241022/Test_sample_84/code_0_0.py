def apply_rules(s):
    # Returns (new_string, changed)
    
    # Rule 1: Replace first 'abc' with 'cab'
    if 'abc' in s:
        pos = s.find('abc')
        return s[:pos] + 'cab' + s[pos+3:], True
    
    # Rule 2: If starts with 'ac', replace with 'zz'
    if s.startswith('ac'):
        return 'zz' + s[2:], True
    
    # Rule 3: Replace first 'acb' with 'bca'
    if 'acb' in s:
        pos = s.find('acb')
        return s[:pos] + 'bca' + s[pos+3:], True
    
    # Rule 4: If even number of 'b's (and at least one), append 'ab'
    b_count = s.count('b')
    if b_count > 0 and b_count % 2 == 0:
        return s + 'ab', True
    
    # Rule 5: If ends with 'cc', replace with 'b' and prepend 'a'
    if s.endswith('cc'):
        return 'a' + s[:-2] + 'b', True
    
    # Rule 6: If ends with 'ba', replace with 'ab'
    if s.endswith('ba'):
        return s[:-2] + 'ab', True
    
    return s, False

# Initial string
s = "abbbaaabcacbcbbcaba"
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