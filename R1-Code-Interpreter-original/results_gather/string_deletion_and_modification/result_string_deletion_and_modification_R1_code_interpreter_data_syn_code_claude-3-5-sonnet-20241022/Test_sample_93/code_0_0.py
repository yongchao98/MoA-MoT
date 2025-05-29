def apply_rules(s):
    # Returns (new_string, changed)
    
    # Rule 1: abc -> cab
    if 'abc' in s:
        pos = s.find('abc')
        return s[:pos] + 'cab' + s[pos+3:], True
    
    # Rule 2: bca -> ''
    if 'bca' in s:
        pos = s.find('bca')
        return s[:pos] + s[pos+3:], True
    
    # Rule 3: acb -> bca
    if 'acb' in s:
        pos = s.find('acb')
        return s[:pos] + 'bca' + s[pos+3:], True
    
    # Rule 4: prefix ca -> bb + c at end
    if s.startswith('ca'):
        return 'bb' + s[2:] + 'c', True
    
    # Rule 5: suffix cc -> a + b
    if s.endswith('cc'):
        return 'a' + s[:-2] + 'b', True
    
    # Rule 6: prefix ab -> ca
    if s.startswith('ab'):
        return 'ca' + s[2:], True
    
    return s, False

def process_string(s):
    steps = [s]
    while True:
        new_s, changed = apply_rules(s)
        if not changed:
            break
        s = new_s
        steps.append(s)
    
    print("Steps:")
    for i, step in enumerate(steps):
        print(f"Step {i}: {step}")
    return s

# Initial string
initial = "cbcbabcabcacbba"
final = process_string(initial)
print(f"\nFinal result: {final}")