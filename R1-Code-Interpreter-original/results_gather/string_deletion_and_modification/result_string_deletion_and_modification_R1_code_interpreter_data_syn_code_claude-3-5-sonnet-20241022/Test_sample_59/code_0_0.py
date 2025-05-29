def apply_rules(s):
    def rule1(s):  # Replace first 'abc' with 'cab'
        if 'abc' in s:
            i = s.index('abc')
            return s[:i] + 'cab' + s[i+3:], True
        return s, False

    def rule2(s):  # Remove first 'ca' after first char
        if len(s) > 1 and 'ca' in s[1:]:
            i = s[1:].index('ca') + 1
            return s[:i] + s[i+2:], True
        return s, False

    def rule3(s):  # Replace ending 'ba' with 'ab'
        if s.endswith('ba'):
            return s[:-2] + 'ab', True
        return s, False

    def rule4(s):  # Remove last char if ends with 'ca'
        if s.endswith('ca'):
            return s[:-1], True
        return s, False

    def rule5(s):  # Append 'ab' if even number of b's
        b_count = s.count('b')
        if b_count > 0 and b_count % 2 == 0:
            return s + 'ab', True
        return s, False

    def rule6(s):  # Remove first char if starts with 'aa'
        if s.startswith('aa'):
            return s[1:], True
        return s, False

    current = s
    steps = [current]
    
    while True:
        changed = False
        
        # Try each rule in order
        for rule in [rule1, rule2, rule3, rule4, rule5, rule6]:
            new_s, was_changed = rule(current)
            if was_changed:
                current = new_s
                steps.append(f"{current}")
                changed = True
                break
        
        if not changed:
            break
    
    print("\n".join(steps))
    return current

initial = "abcbabaccacaba"
result = apply_rules(initial)