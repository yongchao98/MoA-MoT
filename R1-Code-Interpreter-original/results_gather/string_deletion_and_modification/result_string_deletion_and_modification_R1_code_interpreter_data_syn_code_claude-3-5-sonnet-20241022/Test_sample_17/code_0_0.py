def apply_rules(s):
    # Returns (new_string, changed)
    # Rule 1: If suffix is 'ac', replace with 'cb'
    if s.endswith('ac'):
        return (s[:-2] + 'cb', True)
    
    # Rule 2: If prefix is 'bc', delete first two chars and append 'aa'
    if s.startswith('bc'):
        return (s[2:] + 'aa', True)
    
    # Rule 3: If ends with 'ca', remove last character
    if s.endswith('ca'):
        return (s[:-1], True)
    
    # Rule 4: If prefix is 'aa', remove first character
    if s.startswith('aa'):
        return (s[1:], True)
    
    # Rule 5: If starts with 'ac', replace first two chars with 'zz'
    if s.startswith('ac'):
        return ('zz' + s[2:], True)
    
    # Rule 6: If contains 'bca', delete first occurrence
    if 'bca' in s:
        idx = s.index('bca')
        return (s[:idx] + s[idx+3:], True)
    
    return (s, False)

def process_string(initial):
    current = initial
    while True:
        new_string, changed = apply_rules(current)
        if not changed:
            break
        print(f"{current} -> {new_string}")
        current = new_string
    return current

initial = "bcacccbbcabcacc"
final = process_string(initial)
print(f"\nFinal result: {final}")