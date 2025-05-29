def apply_rules(s):
    # Returns (new_string, changed)
    # Rule 1: If ends with 'aa', replace with 'cc'
    if s.endswith('aa'):
        return s[:-2] + 'cc', True
    
    # Rule 2: If ends with 'bb', delete last two characters
    if s.endswith('bb'):
        return s[:-2], True
    
    # Rule 3: If starts with 'cb', replace with 'aa' and delete last char
    if s.startswith('cb'):
        return 'aa' + s[2:-1], True
    
    # Rule 4: If length > 15, remove middle character
    if len(s) > 15:
        mid = len(s) // 2
        return s[:mid] + s[mid+1:], True
    
    # Rule 5: If starts with 'bb', remove second character
    if s.startswith('bb'):
        return s[0] + s[2:], True
    
    # Rule 6: If ends with 'ac', replace with 'cb'
    if s.endswith('ac'):
        return s[:-2] + 'cb', True
    
    return s, False

def process_string(initial):
    current = initial
    steps = [current]
    
    while True:
        new_string, changed = apply_rules(current)
        if not changed:
            break
        current = new_string
        steps.append(current)
    
    print("Initial string:", initial)
    for i, s in enumerate(steps):
        print(f"Step {i}: {s}")
    print("\nFinal string:", current)

# Process the given string
process_string("aaccaaacbbcbaabcbbc")