def apply_rules(s):
    def rule1(s): # ends with 'aa' -> 'cc'
        return s[:-2] + 'cc' if s.endswith('aa') else s
    
    def rule2(s): # ends with 'bb' -> delete last 2
        return s[:-2] if s.endswith('bb') else s
    
    def rule3(s): # starts with 'ab' -> 'ca'
        return 'ca' + s[2:] if s.startswith('ab') else s
    
    def rule4(s): # starts with 'bc' -> delete first 2 and add 'aa'
        return s[2:] + 'aa' if s.startswith('bc') else s
    
    def rule5(s): # contains 'abc' -> replace first with 'cab'
        return s.replace('abc', 'cab', 1) if 'abc' in s else s
    
    def rule6(s): # contains 'bca' -> delete first occurrence
        return s.replace('bca', '', 1) if 'bca' in s else s
    
    current = s
    while True:
        changed = False
        # Try each rule in order
        for rule in [rule1, rule2, rule3, rule4, rule5, rule6]:
            new_s = rule(current)
            if new_s != current:
                print(f"Applied rule, changed {current} to {new_s}")
                current = new_s
                changed = True
                break
        if not changed:
            break
    return current

initial = "cbababccac"
result = apply_rules(initial)
print(f"\nFinal result: {result}")