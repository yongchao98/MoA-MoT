def apply_rules(s):
    # Returns (new_string, changed)
    # Rule 1: prefix 'bc' -> delete first two and append 'aa'
    if s.startswith('bc'):
        return s[2:] + 'aa', True
    
    # Rule 2: prefix 'aa' -> remove first character
    if s.startswith('aa'):
        return s[1:], True
    
    # Rule 3: prefix 'ac' -> replace with 'zz'
    if s.startswith('ac'):
        return 'zz' + s[2:], True
    
    # Rule 4: even number of 'b's -> append 'ab'
    if 'b' in s and s.count('b') % 2 == 0:
        return s + 'ab', True
    
    # Rule 5: prefix 'ab' -> replace with 'ca'
    if s.startswith('ab'):
        return 'ca' + s[2:], True
    
    # Rule 6: prefix 'cb' -> replace with 'aa' and delete last char
    if s.startswith('cb'):
        return 'aa' + s[2:-1], True
    
    return s, False

# Initial string
s = "cbbcbbaabbcb"
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