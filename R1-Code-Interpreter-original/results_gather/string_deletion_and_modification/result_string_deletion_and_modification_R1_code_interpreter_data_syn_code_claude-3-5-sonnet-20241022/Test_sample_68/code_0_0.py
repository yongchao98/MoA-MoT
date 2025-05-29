def apply_rules(s):
    def rule1(s):  # ends with 'ba' -> 'ab'
        if s.endswith('ba'):
            return s[:-2] + 'ab'
        return s

    def rule2(s):  # starts with 'ca' -> 'bb' + 'c'
        if s.startswith('ca'):
            return 'bb' + s[2:] + 'c'
        return s

    def rule3(s):  # starts with 'ab' -> 'ca'
        if s.startswith('ab'):
            return 'ca' + s[2:]
        return s

    def rule4(s):  # starts with 'bb' -> remove second char
        if s.startswith('bb'):
            return 'b' + s[2:]
        return s

    def rule5(s):  # ends with 'bb' -> delete last two chars
        if s.endswith('bb'):
            return s[:-2]
        return s

    def rule6(s):  # contains 'ca' not at start -> remove first occurrence after first char
        if len(s) > 1:
            pos = s[1:].find('ca')
            if pos != -1:
                pos += 1  # adjust for the slice
                return s[:pos] + s[pos+2:]
        return s

    rules = [rule1, rule2, rule3, rule4, rule5, rule6]
    steps = []
    current = s
    steps.append(current)
    
    while True:
        changed = False
        for rule in rules:
            new_s = rule(current)
            if new_s != current:
                current = new_s
                steps.append(current)
                changed = True
                break
        if not changed:
            break
    
    print("Steps:")
    for i, step in enumerate(steps):
        print(f"Step {i}: {step}")
    return current

initial = "abaacaabcab"
final = apply_rules(initial)
print(f"\nFinal result: {final}")