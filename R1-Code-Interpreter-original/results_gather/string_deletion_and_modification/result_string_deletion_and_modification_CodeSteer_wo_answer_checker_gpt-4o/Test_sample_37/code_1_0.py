def transform_string(s):
    while True:
        original = s
        
        # Step 1: If the string starts with 'ac', replace the first two characters with 'zz'.
        if s.startswith('ac'):
            s = 'zz' + s[2:]
        
        # Step 2: If the string ends with 'bb', delete the last two characters.
        elif s.endswith('bb'):
            s = s[:-2]
        
        # Step 3: If the string contains an even number of 'b's (and at least one 'b'), append 'ab' at the end.
        elif s.count('b') % 2 == 0 and s.count('b') > 0:
            s += 'ab'
        
        # Step 4: If the string starts with 'cb', replace it with 'aa' and delete the last character.
        elif s.startswith('cb'):
            s = 'aa' + s[2:-1]
        
        # Step 5: If the string contains 'abc', replace the first occurrence with 'cab'.
        elif 'abc' in s:
            s = s.replace('abc', 'cab', 1)
        
        # Step 6: If the string ends with 'aa', replace it with 'cc'.
        elif s.endswith('aa'):
            s = s[:-2] + 'cc'
        
        # If no changes were made, break the loop
        if s == original:
            break
    
    return s

final_string = transform_string("accabcabccbcccccabca")
print(f'<<< "{final_string}" >>>')