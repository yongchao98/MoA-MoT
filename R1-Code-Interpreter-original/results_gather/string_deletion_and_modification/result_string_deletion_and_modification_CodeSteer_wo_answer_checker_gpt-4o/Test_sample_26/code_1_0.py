def transform_string(s):
    while True:
        original = s
        
        # 1. If the string starts with 'cc', remove the first two characters.
        if s.startswith('cc'):
            s = s[2:]
        
        # 2. If the string contains 'acb', replace the first occurrence with 'bca'.
        elif 'acb' in s:
            s = s.replace('acb', 'bca', 1)
        
        # 3. If the string prefix is 'aa', remove the first character.
        elif s.startswith('aa'):
            s = s[1:]
        
        # 4. If the string ends with 'ca', remove the last character.
        elif s.endswith('ca'):
            s = s[:-1]
        
        # 5. If the string suffix is 'cc', replace it with 'b' and prepend 'a' to the start.
        elif s.endswith('cc'):
            s = 'a' + s[:-2] + 'b'
        
        # 6. If the string contains an even number of 'b's (and at least one 'b'), append 'ab' at the end.
        elif s.count('b') % 2 == 0 and s.count('b') > 0:
            s += 'ab'
        
        # If no changes were made, break the loop
        if s == original:
            break
    
    return s

final_string = transform_string("aaabcbbccbcaacc")
print(f'<<< "{final_string}" >>>')