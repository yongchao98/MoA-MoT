def apply_rules(s):
    def rule1(s):  # ends with 'ba' -> 'ab'
        if s.endswith('ba'):
            return s[:-2] + 'ab', True
        return s, False

    def rule2(s):  # suffix 'ac' -> 'cb'
        if s.endswith('ac'):
            return s[:-2] + 'cb', True
        return s, False

    def rule3(s):  # even number of b's -> append 'ab'
        b_count = s.count('b')
        if b_count > 0 and b_count % 2 == 0:
            return s + 'ab', True
        return s, False

    def rule4(s):  # starts with 'bb' -> remove second char
        if s.startswith('bb'):
            return s[0] + s[2:], True
        return s, False

    def rule5(s):  # length > 15 -> remove middle char
        if len(s) > 15:
            mid = len(s) // 2
            return s[:mid] + s[mid+1:], True
        return s, False

    def rule6(s):  # contains 'abc' -> replace first with 'cab'
        if 'abc' in s:
            idx = s.index('abc')
            return s[:idx] + 'cab' + s[idx+3:], True
        return s, False

    original = s
    steps = [original]
    
    while True:
        changed = False
        
        # Try each rule in order
        for rule in [rule1, rule2, rule3, rule4, rule5, rule6]:
            s, was_changed = rule(s)
            if was_changed:
                steps.append(s)
                changed = True
                break
        
        if not changed:
            break
    
    print("Steps:")
    for i, step in enumerate(steps):
        print(f"Step {i}: {step}")
    
    return s

initial = "caccbacbbaacacaccac"
final = apply_rules(initial)
print(f"\nFinal result: {final}")