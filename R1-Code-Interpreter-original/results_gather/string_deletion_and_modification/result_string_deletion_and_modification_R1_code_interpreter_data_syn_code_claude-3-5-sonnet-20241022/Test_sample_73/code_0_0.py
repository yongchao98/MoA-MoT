def apply_rules(s):
    # Returns (new_string, changed)
    
    # Rule 1: If suffix is 'ac', replace with 'cb'
    if s.endswith('ac'):
        return s[:-2] + 'cb', True
    
    # Rule 2: If prefix is 'aa', remove first character
    if s.startswith('aa'):
        return s[1:], True
    
    # Rule 3: If contains 'acb', replace first occurrence with 'bca'
    if 'acb' in s:
        pos = s.find('acb')
        return s[:pos] + 'bca' + s[pos+3:], True
    
    # Rule 4: If length > 15, remove middle character
    if len(s) > 15:
        mid = len(s) // 2
        return s[:mid] + s[mid+1:], True
    
    # Rule 5: If prefix is 'ab', replace with 'ca'
    if s.startswith('ab'):
        return 'ca' + s[2:], True
    
    # Rule 6: If contains 'ca' (not at start), find first occurrence after first char
    if len(s) > 1 and 'ca' in s[1:]:
        pos = s.find('ca', 1)
        return s[:pos] + s[pos+2:], True
    
    return s, False

def process_string(s):
    steps = 1
    while True:
        print(f"Step {steps}: {s}")
        new_s, changed = apply_rules(s)
        if not changed:
            break
        s = new_s
        steps += 1
    return s

initial = "ccacbcbbcacaaabbabb"
final = process_string(initial)
print(f"\nFinal result: {final}")