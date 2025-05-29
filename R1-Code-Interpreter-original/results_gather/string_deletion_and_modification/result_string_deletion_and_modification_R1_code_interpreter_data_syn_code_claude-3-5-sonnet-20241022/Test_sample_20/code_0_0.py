def apply_rules(s):
    def rule1(s):  # Delete first 'bca'
        if 'bca' in s:
            i = s.index('bca')
            return s[:i] + s[i+3:], True
        return s, False

    def rule2(s):  # If even number of b's, append 'ab'
        if s.count('b') > 0 and s.count('b') % 2 == 0:
            return s + 'ab', True
        return s, False

    def rule3(s):  # If starts with 'aa', remove first char
        if s.startswith('aa'):
            return s[1:], True
        return s, False

    def rule4(s):  # If ends with 'aa', replace with 'cc'
        if s.endswith('aa'):
            return s[:-2] + 'cc', True
        return s, False

    def rule5(s):  # If starts with 'bb', remove second char
        if s.startswith('bb'):
            return s[0] + s[2:], True
        return s, False

    def rule6(s):  # Replace first 'acb' with 'bca'
        if 'acb' in s:
            i = s.index('acb')
            return s[:i] + 'bca' + s[i+3:], True
        return s, False

    current = s
    steps = [current]
    while True:
        changed = False
        
        # Try each rule in order
        for rule in [rule1, rule2, rule3, rule4, rule5, rule6]:
            new_s, rule_changed = rule(current)
            if rule_changed:
                current = new_s
                steps.append(f"{current}")
                changed = True
                break
        
        if not changed:
            break

    print("\n".join(steps))
    return current

initial = "baacbabbabbbc"
final = apply_rules(initial)
print("\nFinal answer:")
print(f'<<< "{final}" >>>')