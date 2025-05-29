def transform_string(s):
    print(f"Start: {s}")
    
    # Store each transformation
    transformations = []
    transformations.append(s)
    
    while True:
        changed = False
        
        # Rule 1: If starts with 'bb', remove second character
        if s.startswith('bb'):
            s = s[0] + s[2:]
            print(f"Rule 1: {s}")
            transformations.append(s)
            changed = True
            continue
        
        # Rule 2: If contains 'acb', replace with 'bca'
        if 'acb' in s:
            i = s.find('acb')
            s = s[:i] + 'bca' + s[i+3:]
            print(f"Rule 2: {s}")
            transformations.append(s)
            changed = True
            continue
        
        # Rule 3: If ends with 'bb', delete last two characters
        if s.endswith('bb'):
            s = s[:-2]
            print(f"Rule 3: {s}")
            transformations.append(s)
            changed = True
            continue
        
        # Rule 4: If ends with 'cc', replace with 'b' and prepend 'a'
        if s.endswith('cc'):
            s = 'a' + s[:-2] + 'b'
            print(f"Rule 4: {s}")
            transformations.append(s)
            changed = True
            continue
        
        # Rule 5: If contains 'abc', replace with 'cab'
        if 'abc' in s:
            i = s.find('abc')
            s = s[:i] + 'cab' + s[i+3:]
            print(f"Rule 5: {s}")
            transformations.append(s)
            changed = True
            continue
        
        # Rule 6: If ends with 'ba', replace with 'ab'
        if s[-2:] == 'ba':  # More explicit check for 'ba' ending
            s = s[:-2] + 'ab'
            print(f"Rule 6: {s}")
            transformations.append(s)
            changed = True
            continue
        
        if not changed:
            break
    
    print("\nAll transformations:")
    for i, t in enumerate(transformations):
        print(f"Step {i}: {t}")
    
    return s

initial = "cbbcabaccabbcbab"
result = transform_string(initial)