def apply_rules(s):
    def debug_print(rule_num, old_s, new_s):
        print(f"Rule {rule_num} applied: {old_s} -> {new_s}")
        return new_s

    def find_abc(s):
        # Look for overlapping matches of 'abc'
        for i in range(len(s)-2):
            if s[i:i+3] == 'abc':
                return i
        return -1

    iteration = 0
    while True:
        print(f"\nIteration {iteration}. Current string: {s}")
        original = s
        
        # Rule 1: ends with 'ba' -> 'ab'
        if s.endswith('ba'):
            s = debug_print(1, s, s[:-2] + 'ab')
            continue
            
        # Rule 2: contains 'bca' -> delete first occurrence
        if 'bca' in s:
            idx = s.find('bca')
            s = debug_print(2, s, s[:idx] + s[idx+3:])
            continue
            
        # Rule 3: starts with 'bb' -> remove second character
        if s.startswith('bb'):
            s = debug_print(3, s, s[0] + s[2:])
            continue
            
        # Rule 4: starts with 'aa' -> remove first character
        if s.startswith('aa'):
            s = debug_print(4, s, s[1:])
            continue
            
        # Rule 5: starts with 'cb' -> replace with 'aa' and delete last char
        if s.startswith('cb'):
            s = debug_print(5, s, 'aa' + s[2:-1])
            continue
            
        # Rule 6: contains 'abc' -> replace first occurrence with 'cab'
        abc_pos = find_abc(s)
        if abc_pos != -1:
            new_s = s[:abc_pos] + 'cab' + s[abc_pos+3:]
            s = debug_print(6, s, new_s)
            continue
            
        if s == original:
            break
            
        iteration += 1

    return s

initial = "aaaaacabacbcbc"
print("Initial string:", initial)
final = apply_rules(initial)
print("\nFinal string:", final)