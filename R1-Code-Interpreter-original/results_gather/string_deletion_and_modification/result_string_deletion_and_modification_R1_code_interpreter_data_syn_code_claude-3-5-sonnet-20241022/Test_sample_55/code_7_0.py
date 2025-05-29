def transform_string(s):
    print(f"Start: {s}")
    
    while True:
        current = s
        print(f"\nChecking string: {s}")
        print(f"Last two characters: '{s[-2:]}'")
        
        # Rule 1: If starts with 'bb', remove second character
        if len(s) >= 2 and s[0:2] == 'bb':
            s = s[0] + s[2:]
            print(f"Applied Rule 1: {s}")
            continue
        
        # Rule 2: If contains 'acb', replace with 'bca'
        if 'acb' in s:
            i = s.find('acb')
            s = s[:i] + 'bca' + s[i+3:]
            print(f"Applied Rule 2: {s}")
            continue
        
        # Rule 3: If ends with 'bb', delete last two characters
        if len(s) >= 2 and s[-2:] == 'bb':
            s = s[:-2]
            print(f"Applied Rule 3: {s}")
            continue
        
        # Rule 4: If ends with 'cc', replace with 'b' and prepend 'a'
        if len(s) >= 2 and s[-2:] == 'cc':
            s = 'a' + s[:-2] + 'b'
            print(f"Applied Rule 4: {s}")
            continue
        
        # Rule 5: If contains 'abc', replace with 'cab'
        if 'abc' in s:
            i = s.find('abc')
            s = s[:i] + 'cab' + s[i+3:]
            print(f"Applied Rule 5: {s}")
            continue
        
        # Rule 6: If ends with 'ba', replace with 'ab'
        if len(s) >= 2 and s[-2:] == 'ba':
            s = s[:-2] + 'ab'
            print(f"Applied Rule 6: {s}")
            continue
        
        # If no rules applied, break
        if s == current:
            print("No rules applied - stopping")
            break
    
    print(f"\nFinal string: {s}")
    return s

s = "cbbcabaccabbcbab"
result = transform_string(s)