def apply_operations(s):
    def try_operations(s):
        # Operation 1: If prefix is 'ab', replace with 'ca'
        if s.startswith('ab'):
            return 'ca' + s[2:]
        
        # Operation 2: If prefix is 'cb', replace with 'aa' and delete last char
        if s.startswith('cb'):
            return 'aa' + s[2:-1]
        
        # Operation 3: If suffix is 'ac', replace with 'cb'
        if s.endswith('ac'):
            return s[:-2] + 'cb'
        
        # Operation 4: If prefix is 'bc', delete first two chars and append 'aa'
        if s.startswith('bc'):
            return s[2:] + 'aa'
        
        # Operation 5: If ends with 'aa', replace with 'cc'
        if s.endswith('aa'):
            return s[:-2] + 'cc'
        
        # Operation 6: If starts with 'bb', remove second character
        if s.startswith('bb'):
            return s[0] + s[2:]
        
        return s

    seen = set()
    current = s
    while True:
        seen.add(current)
        next_s = try_operations(current)
        print(f"Current: {current} -> Next: {next_s}")
        if next_s == current:
            break
        if next_s in seen:
            break
        current = next_s
    
    return current

initial = "bbcaaccbbccbbbaa"
result = apply_operations(initial)
print(f"\nFinal result: {result}")