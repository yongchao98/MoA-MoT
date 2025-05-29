def transform_string(s):
    while True:
        original = s
        
        # Operation 1: If the string suffix is 'ac', replace it with 'cb'.
        if s.endswith('ac'):
            s = s[:-2] + 'cb'
        
        # Operation 2: If the string prefix is 'bc', delete the first two characters and append 'aa' to the end.
        if s.startswith('bc'):
            s = s[2:] + 'aa'
        
        # Operation 3: If the string ends with 'ca', remove the last character.
        if s.endswith('ca'):
            s = s[:-1]
        
        # Operation 4: If the string prefix is 'aa', remove the first character.
        if s.startswith('aa'):
            s = s[1:]
        
        # Operation 5: If the string starts with 'ac', replace the first two characters with 'zz'.
        if s.startswith('ac'):
            s = 'zz' + s[2:]
        
        # Operation 6: If the string contains 'bca', delete the first occurrence entirely.
        if 'bca' in s:
            s = s.replace('bca', '', 1)
        
        # If no changes were made, break the loop
        if s == original:
            break
    
    return s

# Initial string
initial_string = "bcacccbbcabcacc"
final_string = transform_string(initial_string)
print(f'<<< "{final_string}" >>>')