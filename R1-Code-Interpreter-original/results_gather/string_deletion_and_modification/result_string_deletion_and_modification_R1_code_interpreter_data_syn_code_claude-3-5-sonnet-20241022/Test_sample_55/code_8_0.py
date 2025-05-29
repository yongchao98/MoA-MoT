def transform_string(s):
    print(f"Initial: {s}")
    
    def check_and_apply_rules(s):
        # Rule 1: If starts with 'bb', remove second character
        if len(s) >= 2 and s[0:2] == 'bb':
            return s[0] + s[2:], 1
        
        # Rule 2: If contains 'acb', replace with 'bca'
        if 'acb' in s:
            i = s.find('acb')
            return s[:i] + 'bca' + s[i+3:], 2
        
        # Rule 3: If ends with 'bb', delete last two characters
        if len(s) >= 2 and s[-2:] == 'bb':
            return s[:-2], 3
        
        # Rule 4: If ends with 'cc', replace with 'b' and prepend 'a'
        if len(s) >= 2 and s[-2:] == 'cc':
            return 'a' + s[:-2] + 'b', 4
        
        # Rule 5: If contains 'abc', replace with 'cab'
        if 'abc' in s:
            i = s.find('abc')
            return s[:i] + 'cab' + s[i+3:], 5
        
        # Rule 6: If ends with 'ba', replace with 'ab'
        if len(s) >= 2 and s[-2:] == 'ba':
            return s[:-2] + 'ab', 6
        
        return s, 0

    step = 0
    while True:
        new_s, rule = check_and_apply_rules(s)
        if rule > 0:
            step += 1
            print(f"Step {step}: Applied Rule {rule}: {s} -> {new_s}")
            s = new_s
        else:
            break
    
    print(f"Final: {s}")
    return s

# Test with the original string
s = "cbbcabaccabbcbab"
result = transform_string(s)