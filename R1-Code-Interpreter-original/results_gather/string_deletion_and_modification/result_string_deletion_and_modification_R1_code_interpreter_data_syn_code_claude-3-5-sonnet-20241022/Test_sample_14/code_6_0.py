def apply_rules(s):
    def debug_print(rule_num, old_s, new_s, pattern=""):
        print(f"Rule {rule_num} applied: {old_s} -> {new_s}")
        if pattern:
            print(f"Found: {pattern}")
        return new_s

    iteration = 0
    while True:
        print(f"\nIteration {iteration}: {s}")
        
        original = s
        changed = False

        # Rule 1: ends with 'ba' -> 'ab'
        if s.endswith('ba'):
            s = debug_print(1, s, s[:-2] + 'ab', "'ba' at end")
            changed = True
            
        # Rule 2: contains 'bca' -> delete first occurrence
        elif 'bca' in s:
            idx = s.find('bca')
            s = debug_print(2, s, s[:idx] + s[idx+3:], "'bca' found")
            changed = True
            
        # Rule 3: starts with 'bb' -> remove second character
        elif s.startswith('bb'):
            s = debug_print(3, s, s[0] + s[2:], "'bb' at start")
            changed = True
            
        # Rule 4: starts with 'aa' -> remove first character
        elif s.startswith('aa'):
            s = debug_print(4, s, s[1:], "'aa' at start")
            changed = True
            
        # Rule 5: starts with 'cb' -> replace with 'aa' and delete last char
        elif s.startswith('cb'):
            s = debug_print(5, s, 'aa' + s[2:-1], "'cb' at start")
            changed = True
            
        # Rule 6: contains 'abc' -> replace first occurrence with 'cab'
        else:
            idx = s.find('abc')
            if idx != -1:
                s = debug_print(6, s, s[:idx] + 'cab' + s[idx+3:], "'abc' found")
                changed = True
        
        if not changed:
            break
            
        iteration += 1

    return s

initial = "aaaaacabacbcbc"
print("Initial string:", initial)
final = apply_rules(initial)
print("\nFinal string:", final)