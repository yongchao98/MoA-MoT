def apply_operations(s):
    def op1(s):  # If starts with 'bb', remove second character
        if s.startswith('bb'):
            return s[0] + s[2:], True
        return s, False
    
    def op2(s):  # If contains 'acb', replace first occurrence with 'bca'
        if 'acb' in s:
            idx = s.index('acb')
            return s[:idx] + 'bca' + s[idx+3:], True
        return s, False
    
    def op3(s):  # If suffix is 'bb', delete last two characters
        if s.endswith('bb'):
            return s[:-2], True
        return s, False
    
    def op4(s):  # If suffix is 'cc', replace with 'b' and prepend 'a'
        if s.endswith('cc'):
            return 'a' + s[:-2] + 'b', True
        return s, False
    
    def op5(s):  # If contains 'abc', replace first occurrence with 'cab'
        if 'abc' in s:
            idx = s.index('abc')
            return s[:idx] + 'cab' + s[idx+3:], True
        return s, False
    
    def op6(s):  # If ends with 'ba', replace with 'ab'
        if len(s) >= 2 and s[-2:] == 'ba':
            return s[:-2] + 'ab', True
        return s, False

    current = s
    print(f"Initial string: {current}")
    
    while True:
        changed = False
        # Try each operation
        if current.startswith('bb'):
            current = current[0] + current[2:]
            print(f"Rule 1: {current}")
            changed = True
            continue
            
        if 'acb' in current:
            idx = current.index('acb')
            current = current[:idx] + 'bca' + current[idx+3:]
            print(f"Rule 2: {current}")
            changed = True
            continue
            
        if current.endswith('bb'):
            current = current[:-2]
            print(f"Rule 3: {current}")
            changed = True
            continue
            
        if current.endswith('cc'):
            current = 'a' + current[:-2] + 'b'
            print(f"Rule 4: {current}")
            changed = True
            continue
            
        if 'abc' in current:
            idx = current.index('abc')
            current = current[:idx] + 'cab' + current[idx+3:]
            print(f"Rule 5: {current}")
            changed = True
            continue
            
        if current.endswith('ba'):
            current = current[:-2] + 'ab'
            print(f"Rule 6: {current}")
            changed = True
            continue
            
        if not changed:
            break
    
    print(f"Final string: {current}")
    return current

result = apply_operations("cbbcabaccabbcbab")