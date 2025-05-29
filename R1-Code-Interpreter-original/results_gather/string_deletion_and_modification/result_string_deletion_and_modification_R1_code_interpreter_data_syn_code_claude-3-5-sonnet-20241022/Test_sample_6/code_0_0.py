def apply_rules(s):
    # Returns (new_string, changed)
    
    # Rule 1: If prefix is 'bc', delete first two and append 'aa'
    if s.startswith('bc'):
        return s[2:] + 'aa', True
    
    # Rule 2: If contains 'abc', replace first with 'cab'
    if 'abc' in s:
        pos = s.find('abc')
        return s[:pos] + 'cab' + s[pos+3:], True
    
    # Rule 3: If suffix is 'ac', replace with 'cb'
    if s.endswith('ac'):
        return s[:-2] + 'cb', True
    
    # Rule 4: If prefix is 'ab', replace with 'ca'
    if s.startswith('ab'):
        return 'ca' + s[2:], True
    
    # Rule 5: If contains 'acb', replace with 'bca'
    if 'acb' in s:
        pos = s.find('acb')
        return s[:pos] + 'bca' + s[pos+3:], True
    
    # Rule 6: If suffix is 'bb', delete last two characters
    if s.endswith('bb'):
        return s[:-2], True
    
    return s, False

# Initial string
s = "acbbaacabcacbbccaaa"

# Keep track of all states to detect cycles
steps = [s]
while True:
    s, changed = apply_rules(s)
    if not changed:
        break
    print(f"After transformation: {s}")

print("\nFinal string:", s)