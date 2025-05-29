def apply_rules(s):
    def find_abc(s):
        # Special check for 'abc' pattern, including when split across other patterns
        for i in range(len(s)-2):
            if s[i] == 'a' and s[i+1] == 'b' and s[i+2] == 'c':
                return i
        return -1

    print(f"Starting with: {s}")
    while True:
        original = s
        changed = False
        
        # Print current state and highlight potential 'abc' patterns
        print(f"\nCurrent string: {s}")
        abc_pos = find_abc(s)
        if abc_pos != -1:
            print(f"Found 'abc' at position {abc_pos}: {s[:abc_pos]}[{s[abc_pos:abc_pos+3]}]{s[abc_pos+3:]}")
        
        # Apply rules in order
        if s.endswith('ba'):
            s = s[:-2] + 'ab'
            print("Rule 1: ends with 'ba' -> 'ab'")
            changed = True
        elif 'bca' in s:
            idx = s.find('bca')
            s = s[:idx] + s[idx+3:]
            print("Rule 2: contains 'bca' -> delete")
            changed = True
        elif s.startswith('bb'):
            s = s[0] + s[2:]
            print("Rule 3: starts with 'bb' -> remove second char")
            changed = True
        elif s.startswith('aa'):
            s = s[1:]
            print("Rule 4: starts with 'aa' -> remove first char")
            changed = True
        elif s.startswith('cb'):
            s = 'aa' + s[2:-1]
            print("Rule 5: starts with 'cb' -> 'aa' and remove last char")
            changed = True
        else:
            abc_pos = find_abc(s)
            if abc_pos != -1:
                s = s[:abc_pos] + 'cab' + s[abc_pos+3:]
                print(f"Rule 6: 'abc' -> 'cab' at position {abc_pos}")
                changed = True
        
        if not changed:
            break
    
    return s

s = "aaaaacabacbcbc"
final = apply_rules(s)
print(f"\nFinal result: {final}")