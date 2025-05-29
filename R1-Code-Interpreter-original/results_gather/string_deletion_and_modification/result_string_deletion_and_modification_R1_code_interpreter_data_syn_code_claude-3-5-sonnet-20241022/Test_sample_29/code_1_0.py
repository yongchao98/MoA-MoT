def apply_rules(s):
    # Returns (new_string, changed, rule_applied)
    # Rule 1: suffix 'ac' -> 'cb'
    if s.endswith('ac'):
        return (s[:-2] + 'cb', True, 1)
    
    # Rule 2: prefix 'bc' -> delete first 2 chars and append 'aa'
    if s.startswith('bc'):
        return (s[2:] + 'aa', True, 2)
    
    # Rule 3: suffix 'bb' -> delete last 2 chars
    if s.endswith('bb'):
        return (s[:-2], True, 3)
    
    # Rule 4: prefix 'ab' -> 'ca'
    if s.startswith('ab'):
        return ('ca' + s[2:], True, 4)
    
    # Rule 5: suffix 'cc' -> prepend 'a' and append 'b'
    if s.endswith('cc'):
        return ('a' + s[:-2] + 'b', True, 5)
    
    # Rule 6: prefix 'cb' -> 'aa' and delete last char
    if s.startswith('cb'):
        return ('aa' + s[2:-1], True, 6)
    
    return (s, False, 0)

def solve(initial_string):
    current = initial_string
    steps = [(current, 0)]
    
    while True:
        new_string, changed, rule = apply_rules(current)
        if not changed:
            break
        current = new_string
        steps.append((current, rule))
    
    print("Steps:")
    for i, (step, rule) in enumerate(steps):
        print(f"Step {i}: {step} (Rule {rule})")
    return current

# Continue from where we left off
initial = "aacabbabcccaacb"
final = solve(initial)
print(f"\nFinal string: {final}")