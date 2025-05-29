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
        
        # Rule 5: If the string starts with 'bc', delete first two chars and append 'aa'
        if s.startswith('bc'):
            return (s[2:] + 'aa', 5)
        
        # Rule 6: If the string ends with 'cc', replace with 'b' and prepend 'a'
        if s.endswith('cc'):
            return ('a' + s[:-2] + 'b', 6)
        
        return (s, 0)  # No rule applied

    current = s
    step = 0
    print(f"Step {step}: {current}")
    
    while True:
        new_s, rule = try_rule(current)
        if new_s == current:  # No rule was applied
            break
            
        step += 1
        print(f"Step {step}: Applied Rule {rule}: {current} -> {new_s}")
        current = new_s
    
    return current

# Let's start with the original string
s = "ccccaaaaccbabcbabcbb"  # Going back to the original string
result = apply_rules(s)
print(f"\nFinal result: {result}")