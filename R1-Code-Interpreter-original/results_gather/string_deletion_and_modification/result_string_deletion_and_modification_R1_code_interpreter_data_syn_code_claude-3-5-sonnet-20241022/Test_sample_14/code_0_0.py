def apply_rules(s):
    def rule1(s):  # ends with 'ba' -> 'ab'
        if s.endswith('ba'):
            return s[:-2] + 'ab', True
        return s, False

    def rule2(s):  # contains 'bca' -> delete first occurrence
        if 'bca' in s:
            idx = s.find('bca')
            return s[:idx] + s[idx+3:], True
        return s, False

    def rule3(s):  # starts with 'bb' -> remove second char
        if s.startswith('bb'):
            return s[0] + s[2:], True
        return s, False

    def rule4(s):  # starts with 'aa' -> remove first char
        if s.startswith('aa'):
            return s[1:], True
        return s, False

    def rule5(s):  # starts with 'cb' -> replace with 'aa' and delete last char
        if s.startswith('cb'):
            return 'aa' + s[2:-1], True
        return s, False

    def rule6(s):  # contains 'abc' -> replace first occurrence with 'cab'
        if 'abc' in s:
            idx = s.find('abc')
            return s[:idx] + 'cab' + s[idx+3:], True
        return s, False

    rules = [rule1, rule2, rule3, rule4, rule5, rule6]
    
    current = s
    steps = []
    steps.append(current)
    
    while True:
        any_change = False
        for rule in rules:
            new_s, changed = rule(current)
            if changed:
                current = new_s
                steps.append(current)
                any_change = True
                break
        if not any_change:
            break
    
    print("Steps:")
    for i, step in enumerate(steps):
        print(f"Step {i}: {step}")
    return current

initial = "aaaaacabacbcbc"
final = apply_rules(initial)
print(f"\nFinal result: {final}")