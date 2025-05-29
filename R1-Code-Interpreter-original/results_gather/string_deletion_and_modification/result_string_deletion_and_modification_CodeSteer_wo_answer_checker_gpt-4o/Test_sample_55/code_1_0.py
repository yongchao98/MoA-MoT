def transform_string(s):
    while True:
        original = s
        
        # 1. If the string starts with 'bb', remove the second character.
        if s.startswith('bb'):
            s = s[0] + s[2:]
        
        # 2. If the string contains 'acb', replace the first occurrence with 'bca'.
        if 'acb' in s:
            s = s.replace('acb', 'bca', 1)
        
        # 3. If the string suffix is 'bb', delete the last two characters.
        if s.endswith('bb'):
            s = s[:-2]
        
        # 4. If the string suffix is 'cc', replace it with 'b' and prepend 'a' to the start.
        if s.endswith('cc'):
            s = 'a' + s[:-2] + 'b'
        
        # 5. If the string contains 'abc', replace the first occurrence with 'cab'.
        if 'abc' in s:
            s = s.replace('abc', 'cab', 1)
        
        # 6. If the string ends with 'ba', replace it with 'ab'.
        if s.endswith('ba'):
            s = s[:-2] + 'ab'
        
        # If no changes were made, break the loop
        if s == original:
            break
    
    return s

final_string = transform_string("cbbaacbabccbcbba")
print(f'<<< "{final_string}" >>>')