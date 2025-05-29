def apply_rules(s):
    def check_and_apply(s):
        # Rule 1: ends with 'ba' -> 'ab'
        if s.endswith('ba'):
            return s[:-2] + 'ab', 1
        
        # Rule 2: starts with 'aa' -> remove first 'a'
        if s.startswith('aa'):
            return s[1:], 2
        
        # Rule 3: ends with 'ac' -> 'cb'
        if s.endswith('ac'):
            return s[:-2] + 'cb', 3
        
        # Rule 4: starts with 'cc' -> remove 'cc'
        if s.startswith('cc'):
            return s[2:], 4
        
        # Rule 5: starts with 'bc' -> remove 'bc' and add 'aa'
        if s.startswith('bc'):
            return s[2:] + 'aa', 5
        
        # Rule 6: ends with 'cc' -> add 'a' at start and replace 'cc' with 'b'
        if s.endswith('cc'):
            return 'a' + s[:-2] + 'b', 6
        
        return s, 0

    current = s
    print(f"Start: {current}")
    
    step = 1
    while True:
        new_s, rule = check_and_apply(current)
        if new_s == current:  # No rule applied
            break
            
        print(f"Step {step}: Applied Rule {rule}: {current} -> {new_s}")
        current = new_s
        step += 1
    
    return current

# Test with the string
s = "accbabcbabcbb"
result = apply_rules(s)
print(f"\nFinal result: {result}")