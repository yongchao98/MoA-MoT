def apply_rules(s):
    def print_rule_applied(old_s, new_s, rule_num):
        print(f"Rule {rule_num}: {old_s} -> {new_s}")
        return new_s

    while True:
        original = s
        changed = False
        
        # Rule 1: ends with 'ba' -> 'ab'
        if s.endswith('ba'):
            s = print_rule_applied(s, s[:-2] + 'ab', 1)
            continue

        # Rule 2: contains 'bca' -> delete first occurrence
        if 'bca' in s:
            idx = s.find('bca')
            s = print_rule_applied(s, s[:idx] + s[idx+3:], 2)
            continue

        # Rule 3: starts with 'bb' -> remove second character
        if s.startswith('bb'):
            s = print_rule_applied(s, s[0] + s[2:], 3)
            continue

        # Rule 4: starts with 'aa' -> remove first character
        if s.startswith('aa'):
            s = print_rule_applied(s, s[1:], 4)
            continue

        # Rule 5: starts with 'cb' -> replace with 'aa' and delete last char
        if s.startswith('cb'):
            s = print_rule_applied(s, 'aa' + s[2:-1], 5)
            continue

        # Rule 6: contains 'abc' -> replace first occurrence with 'cab'
        for i in range(len(s)-2):
            if s[i:i+3] == 'abc':
                s = print_rule_applied(s, s[:i] + 'cab' + s[i+3:], 6)
                changed = True
                break
        
        if not changed and s == original:
            break

    return s

s = "aaaaacabacbcbc"
print("Initial:", s)
final = apply_rules(s)
print("Final:", final)