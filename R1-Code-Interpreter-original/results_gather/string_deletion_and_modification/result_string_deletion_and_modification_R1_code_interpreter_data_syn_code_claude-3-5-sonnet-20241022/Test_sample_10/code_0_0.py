def apply_rules(s):
    def rule1(s):  # Replace 'ba' suffix with 'ab'
        if s.endswith('ba'):
            return s[:-2] + 'ab', True
        return s, False

    def rule2(s):  # Remove first char if prefix is 'aa'
        if s.startswith('aa'):
            return s[1:], True
        return s, False

    def rule3(s):  # Replace 'ac' suffix with 'cb'
        if s.endswith('ac'):
            return s[:-2] + 'cb', True
        return s, False

    def rule4(s):  # Remove first two chars if prefix is 'cc'
        if s.startswith('cc'):
            return s[2:], True
        return s, False

    def rule5(s):  # If prefix is 'bc', delete first two chars and append 'aa'
        if s.startswith('bc'):
            return s[2:] + 'aa', True
        return s, False

    def rule6(s):  # If suffix is 'cc', replace with 'b' and prepend 'a'
        if s.endswith('cc'):
            return 'a' + s[:-2] + 'b', True
        return s, False

    current = s
    steps = []
    steps.append(current)
    
    while True:
        changed = False
        
        # Try each rule in order
        for rule_func in [rule1, rule2, rule3, rule4, rule5, rule6]:
            new_s, applied = rule_func(current)
            if applied:
                current = new_s
                steps.append(current)
                changed = True
                break
        
        if not changed:
            break
    
    print("Initial string:", s)
    print("\nSteps:")
    for i, step in enumerate(steps):
        print(f"Step {i}: {step}")
    print("\nFinal string:", current)

# Test with the given string
initial_string = "ccccaaaaccbabcbabcbb"
apply_rules(initial_string)