def count_b(s):
    return s.count('b')

def apply_rules(s):
    # Rule 1: Even number of b's (at least one), append 'ab'
    if count_b(s) > 0 and count_b(s) % 2 == 0:
        return s + 'ab', True, 1
    
    # Rule 2: Prefix 'ab' -> 'ca'
    if s.startswith('ab'):
        return 'ca' + s[2:], True, 2
    
    # Rule 3: Prefix 'cb' -> 'aa' and delete last char
    if s.startswith('cb'):
        return 'aa' + s[2:][:-1], True, 3
    
    # Rule 4: Prefix 'ca' -> 'bb' and append 'c'
    if s.startswith('ca'):
        return 'bb' + s[2:] + 'c', True, 4
    
    # Rule 5: Ends with 'ca' -> remove last char
    if s.endswith('ca'):
        return s[:-1], True, 5
    
    # Rule 6: Contains 'abc' -> replace first occurrence with 'cab'
    if 'abc' in s:
        idx = s.index('abc')
        return s[:idx] + 'cab' + s[idx+3:], True, 6
    
    return s, False, 0

def solve(initial):
    current = initial
    changed = True
    steps = []
    steps.append(current)
    
    while changed:
        changed = False
        current_before = current
        for _ in range(6):  # Try each rule
            current, rule_applied, rule_num = apply_rules(current)
            if rule_applied:
                steps.append(f"{current} (Rule {rule_num})")
                changed = True
                break
    
    print("Steps:")
    for i, step in enumerate(steps):
        print(f"Step {i}: {step}")
    print("\nFinal string:", current)

# Initial string
initial = "bbaccabbcabacabbbbbcab"
solve(initial)