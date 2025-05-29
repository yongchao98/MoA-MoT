def apply_rules(s):
    # Returns (new_string, changed)
    # Rule 1: suffix 'ac' -> 'cb'
    if s.endswith('ac'):
        return (s[:-2] + 'cb', True)
    
    # Rule 2: prefix 'bc' -> delete first 2 chars and append 'aa'
    if s.startswith('bc'):
        return (s[2:] + 'aa', True)
    
    # Rule 3: suffix 'bb' -> delete last 2 chars
    if s.endswith('bb'):
        return (s[:-2], True)
    
    # Rule 4: prefix 'ab' -> 'ca'
    if s.startswith('ab'):
        return ('ca' + s[2:], True)
    
    # Rule 5: suffix 'cc' -> prepend 'a' and append 'b'
    if s.endswith('cc'):
        return ('a' + s[:-2] + 'b', True)
    
    # Rule 6: prefix 'cb' -> 'aa' and delete last char
    if s.startswith('cb'):
        return ('aa' + s[2:-1], True)
    
    return (s, False)

def solve(initial_string):
    current = initial_string
    steps = [current]
    
    while True:
        new_string, changed = apply_rules(current)
        if not changed:
            break
        current = new_string
        steps.append(current)
    
    print("Steps:")
    for i, step in enumerate(steps):
        print(f"Step {i}: {step}")
    return current

# Initial string
initial = "cbcabbabcccaaaacbb"
final = solve(initial)
print(f"\nFinal string: {final}")