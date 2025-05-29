def apply_rules(s):
    # Rule 1: If prefix is 'bc', delete first two and append 'aa'
    if s.startswith('bc'):
        return s[2:] + 'aa', True, "Rule 1"
    
    # Rule 2: If contains 'abc', replace first with 'cab'
    if 'abc' in s:
        pos = s.find('abc')
        return s[:pos] + 'cab' + s[pos+3:], True, "Rule 2"
    
    # Rule 3: If suffix is 'ac', replace with 'cb'
    if s.endswith('ac'):
        return s[:-2] + 'cb', True, "Rule 3"
    
    # Rule 4: If prefix is 'ab', replace with 'ca'
    if s.startswith('ab'):
        pos = s.find('ab')
        return 'ca' + s[2:], True, "Rule 4"
    
    # Rule 5: If contains 'acb', replace with 'bca'
    if 'acb' in s:
        pos = s.find('acb')
        return s[:pos] + 'bca' + s[pos+3:], True, "Rule 5"
    
    # Rule 6: If suffix is 'bb', delete last two characters
    if s.endswith('bb'):
        return s[:-2], True, "Rule 6"
    
    return s, False, "No rule"

# Initial string
s = "caaaccabbcccabaaaaa"
iteration = 0
seen_states = set([s])

while True:
    iteration += 1
    s_old = s
    s, changed, rule = apply_rules(s)
    if changed:
        print(f"Step {iteration}: {rule} -> {s}")
        if s in seen_states:
            print("Cycle detected!")
            break
        seen_states.add(s)
    else:
        break

print("\nFinal string:", s)