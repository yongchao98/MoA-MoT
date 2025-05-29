def apply_rules(s):
    def try_rule(s):
        # Rule 1: If the string ends with 'ba', replace it with 'ab'
        if s.endswith('ba'):
            return (s[:-2] + 'ab', 1)
        
        # Rule 2: If the string starts with 'aa', remove first character
        if s.startswith('aa'):
            return (s[1:], 2)
        
        # Rule 3: If the string ends with 'ac', replace it with 'cb'
        if s.endswith('ac'):
            return (s[:-2] + 'cb', 3)
        
        # Rule 4: If the string starts with 'cc', remove first two characters
        if s.startswith('cc'):
            return (s[2:], 4)
        
        # Rule 5: If the string contains 'bc', process leftmost occurrence
        if 'bc' in s:
            return (s.replace('bc', '', 1) + 'aa', 5)
        
        # Rule 6: If the string ends with 'cc', replace with 'b' and prepend 'a'
        if s.endswith('cc'):
            return ('a' + s[:-2] + 'b', 6)
        
        return (s, 0)

    current = s
    step = 0
    print(f"Step {step}: {current}")
    
    while True:
        old = current
        new_s, rule = try_rule(current)
        if new_s == old:
            break
        step += 1
        print(f"Step {step}: Applied Rule {rule}: {current} -> {new_s}")
        current = new_s
    
    return current

# Test with our string
s = "accbabcbabcbb"
result = apply_rules(s)
print(f"\nFinal result: {result}")