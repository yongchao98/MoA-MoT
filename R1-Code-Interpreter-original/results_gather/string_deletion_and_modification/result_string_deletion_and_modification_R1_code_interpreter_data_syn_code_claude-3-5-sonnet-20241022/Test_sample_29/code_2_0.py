def transform_string(s):
    print(f"Current string: {s}")
    
    # Rule 1: suffix 'ac' -> 'cb'
    if s.endswith('ac'):
        print("Applying Rule 1: suffix 'ac' -> 'cb'")
        return s[:-2] + 'cb'
    
    # Rule 2: prefix 'bc' -> delete first 2 chars and append 'aa'
    if s.startswith('bc'):
        print("Applying Rule 2: prefix 'bc' -> delete and append 'aa'")
        return s[2:] + 'aa'
    
    # Rule 3: suffix 'bb' -> delete last 2 chars
    if s.endswith('bb'):
        print("Applying Rule 3: suffix 'bb' -> delete")
        return s[:-2]
    
    # Rule 4: if 'ab' appears anywhere -> replace with 'ca'
    if 'ab' in s:
        idx = s.index('ab')
        print(f"Applying Rule 4: 'ab' -> 'ca' at position {idx}")
        return s[:idx] + 'ca' + s[idx+2:]
    
    # Rule 5: suffix 'cc' -> prepend 'a' and append 'b'
    if s.endswith('cc'):
        print("Applying Rule 5: suffix 'cc' -> prepend 'a' and append 'b'")
        return 'a' + s[:-2] + 'b'
    
    # Rule 6: prefix 'cb' -> 'aa' and delete last char
    if s.startswith('cb'):
        print("Applying Rule 6: prefix 'cb' -> 'aa' and delete last char")
        return 'aa' + s[2:-1]
    
    return s

def solve(s):
    steps = [s]
    while True:
        new_s = transform_string(s)
        if new_s == s:  # No change occurred
            break
        steps.append(new_s)
        s = new_s
    
    print("\nTransformation steps:")
    for i, step in enumerate(steps):
        print(f"Step {i}: {step}")
    return s

initial = "aacabbabcccaacb"
final = solve(initial)
print(f"\nFinal string: {final}")