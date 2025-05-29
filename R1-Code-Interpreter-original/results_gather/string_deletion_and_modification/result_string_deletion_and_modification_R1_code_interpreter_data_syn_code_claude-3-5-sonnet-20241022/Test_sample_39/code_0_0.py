def apply_rules(s):
    # Returns (new_string, changed)
    
    # Rule 1: If prefix is 'bc', delete first two chars and append 'aa'
    if s.startswith('bc'):
        return s[2:] + 'aa', True
    
    # Rule 2: If contains 'bca', delete first occurrence
    if 'bca' in s:
        i = s.index('bca')
        return s[:i] + s[i+3:], True
    
    # Rule 3: If starts with 'ac', replace first two chars with 'zz'
    if s.startswith('ac'):
        return 'zz' + s[2:], True
    
    # Rule 4: If contains 'abc', replace first occurrence with 'cab'
    if 'abc' in s:
        i = s.index('abc')
        return s[:i] + 'cab' + s[i+3:], True
    
    # Rule 5: If prefix is 'ab', replace with 'ca'
    if s.startswith('ab'):
        return 'ca' + s[2:], True
    
    # Rule 6: If suffix is 'ac', replace with 'cb'
    if s.endswith('ac'):
        return s[:-2] + 'cb', True
    
    return s, False

# Initial string
s = "bcbababbcacbbabaaac"

# Keep track of all transformations
print(f"Start: {s}")
step = 1

while True:
    new_s, changed = apply_rules(s)
    if not changed:
        break
    print(f"Step {step}: {new_s}")
    s = new_s
    step += 1

print(f"Final: {s}")