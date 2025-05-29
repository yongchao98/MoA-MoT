def apply_rules(s):
    # Returns (new_string, changed)
    
    # Rule 1: If suffix is 'cc', replace with 'b' and prepend 'a'
    if s.endswith('cc'):
        return 'a' + s[:-2] + 'b', True
    
    # Rule 2: If prefix is 'ca', replace with 'bb' and append 'c'
    if s.startswith('ca'):
        return 'bb' + s[2:] + 'c', True
    
    # Rule 3: If contains 'abc', replace first occurrence with 'cab'
    if 'abc' in s:
        pos = s.find('abc')
        return s[:pos] + 'cab' + s[pos+3:], True
    
    # Rule 4: If starts with 'bb', remove second character
    if s.startswith('bb'):
        return s[0] + s[2:], True
    
    # Rule 5: If contains 'acb', replace first occurrence with 'bca'
    if 'acb' in s:
        pos = s.find('acb')
        return s[:pos] + 'bca' + s[pos+3:], True
    
    # Rule 6: If contains 'ca' (not at start), look after first character
    if 'ca' in s[1:]:
        pos = s[1:].find('ca') + 1  # Add 1 because we started from index 1
        return s[:pos] + s[pos+2:], True
    
    return s, False

# Initial string
s = "bcbacacabbcaaacbcbc"
steps = 1

while True:
    new_s, changed = apply_rules(s)
    if not changed:
        break
    print(f"Step {steps}: {s} -> {new_s}")
    s = new_s
    steps += 1

print(f"\nFinal string: {s}")