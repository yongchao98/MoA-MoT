def apply_rules(s):
    def rule1(s):  # Replace 'acb' with 'bca'
        if 'acb' in s:
            idx = s.index('acb')
            return s[:idx] + 'bca' + s[idx+3:], True
        return s, False

    def rule2(s):  # Remove middle char if length > 15
        if len(s) > 15:
            mid = len(s) // 2
            return s[:mid] + s[mid+1:], True
        return s, False

    def rule3(s):  # Replace 'abc' with 'cab'
        if 'abc' in s:
            idx = s.index('abc')
            return s[:idx] + 'cab' + s[idx+3:], True
        return s, False

    def rule4(s):  # Append 'ab' if even number of b's
        b_count = s.count('b')
        if b_count > 0 and b_count % 2 == 0:
            return s + 'ab', True
        return s, False

    def rule5(s):  # Delete last two chars if ends with 'bb'
        if s.endswith('bb'):
            return s[:-2], True
        return s, False

    def rule6(s):  # Delete first 'bca'
        if 'bca' in s:
            idx = s.index('bca')
            return s[:idx] + s[idx+3:], True
        return s, False

    current = s
    step = 0
    print(f"Start: {current}")
    
    while True:
        changed = False
        
        # Try each rule in order
        for rule_num, rule in enumerate([rule1, rule2, rule3, rule4, rule5, rule6], 1):
            new_s, was_changed = rule(current)
            if was_changed:
                step += 1
                print(f"Step {step} (Rule {rule_num}): {current} -> {new_s}")
                current = new_s
                changed = True
                break
        
        if not changed:
            break
    
    return current

initial = "abacabaacbccabacaaa"
final = apply_rules(initial)
print(f"\nFinal result: {final}")