def apply_rules(s):
    # Keep track of changes
    original = s
    
    # Rule 1: If ends with 'aa', replace with 'cc'
    if s.endswith('aa'):
        s = s[:-2] + 'cc'
    # Rule 2: If starts with 'cc', remove first two chars
    elif s.startswith('cc'):
        s = s[2:]
    # Rule 3: If contains 'abc', replace first occurrence with 'cab'
    elif 'abc' in s:
        s = s.replace('abc', 'cab', 1)
    # Rule 4: If contains 'bca', delete first occurrence
    elif 'bca' in s:
        s = s.replace('bca', '', 1)
    # Rule 5: If starts with 'ac', replace first two chars with 'zz'
    elif s.startswith('ac'):
        s = 'zz' + s[2:]
    # Rule 6: If contains 'ca' (not at start), find first occurrence after first char
    elif 'ca' in s[1:]:
        idx = s[1:].find('ca') + 1
        s = s[:idx] + s[idx+2:]
    
    return s, s != original

# Initial string
s = "acacbccaabaacbbaa"
step = 1

print(f"Start: {s}")
while True:
    new_s, changed = apply_rules(s)
    if not changed:
        break
    print(f"Step {step}: {new_s}")
    s = new_s
    step += 1

print(f"Final: {s}")