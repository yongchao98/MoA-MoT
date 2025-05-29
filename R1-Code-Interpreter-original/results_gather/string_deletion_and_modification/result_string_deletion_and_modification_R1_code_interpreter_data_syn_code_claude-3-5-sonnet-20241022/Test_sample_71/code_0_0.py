def apply_rules(s):
    # Returns (new_string, changed)
    
    # Rule 1: If prefix is 'ab', replace with 'ca'
    if s.startswith('ab'):
        return 'ca' + s[2:], True
    
    # Rule 2: If starts with 'bb', remove second character
    if s.startswith('bb'):
        return 'b' + s[2:], True
    
    # Rule 3: If contains 'abc', replace first occurrence with 'cab'
    if 'abc' in s:
        pos = s.find('abc')
        return s[:pos] + 'cab' + s[pos+3:], True
    
    # Rule 4: If contains 'bca', delete first occurrence
    if 'bca' in s:
        pos = s.find('bca')
        return s[:pos] + s[pos+3:], True
    
    # Rule 5: If contains 'acb', replace with 'bca'
    if 'acb' in s:
        pos = s.find('acb')
        return s[:pos] + 'bca' + s[pos+3:], True
    
    # Rule 6: If ends with 'aa', replace with 'cc'
    if s.endswith('aa'):
        return s[:-2] + 'cc', True
    
    return s, False

# Initial string
s = "bacccbbcacaaaacb"
steps = [s]

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