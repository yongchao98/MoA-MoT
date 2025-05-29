def apply_rules(s):
    def check_rules(s):
        # Rule 1: If the string ends with 'ba', replace it with 'ab'
        if s.endswith('ba'):
            return s[:-2] + 'ab', True, 1
            
        # Rule 2: If the string prefix is 'aa', remove the first character
        if s.startswith('aa'):
            return s[1:], True, 2
            
        # Rule 3: If the string suffix is 'ac', replace it with 'cb'
        if s.endswith('ac'):
            return s[:-2] + 'cb', True, 3
            
        # Rule 4: If the string starts with 'cc', remove the first two characters
        if s.startswith('cc'):
            return s[2:], True, 4
            
        # Rule 5: If the string contains 'bc' anywhere, process the leftmost occurrence
        if 'bc' in s:
            idx = s.index('bc')
            return s[:idx] + s[idx+2:] + 'aa', True, 5
            
        # Rule 6: If the string suffix is 'cc', replace it with 'b' and prepend 'a'
        if s.endswith('cc'):
            return 'a' + s[:-2] + 'b', True, 6
            
        return s, False, 0

    current = s
    step = 0
    print(f"Step {step}: {current}")
    
    while True:
        new_s, changed, rule = check_rules(current)
        if not changed:
            break
        step += 1
        print(f"Step {step}: Applied Rule {rule}: {current} -> {new_s}")
        current = new_s
    
    return current

# Test with the string
s = "accbabcbabcbb"
result = apply_rules(s)
print(f"\nFinal result: {result}")