def transform_string(s):
    while True:
        original = s
        # 1. If the string ends with 'aa', replace it with 'cc'.
        if s.endswith('aa'):
            s = s[:-2] + 'cc'
        
        # 2. If the string suffix is 'bb', delete the last two characters.
        elif s.endswith('bb'):
            s = s[:-2]
        
        # 3. If the string prefix is 'cb', replace it with 'aa' and delete the last character.
        elif s.startswith('cb'):
            s = 'aa' + s[2:-1]
        
        # 4. If the string length is greater than 15, remove the middle character.
        elif len(s) > 15:
            mid_index = len(s) // 2
            s = s[:mid_index] + s[mid_index+1:]
        
        # 5. If the string starts with 'bb', remove the second character.
        elif s.startswith('bb'):
            s = s[0] + s[2:]
        
        # 6. If the string suffix is 'ac', replace it with 'cb'.
        elif s.endswith('ac'):
            s = s[:-2] + 'cb'
        
        # If no changes were made, break the loop
        if s == original:
            break
    
    return s

final_string = transform_string("aaccaaacbbcbaabcbbc")
print(f'<<< "{final_string}" >>>')