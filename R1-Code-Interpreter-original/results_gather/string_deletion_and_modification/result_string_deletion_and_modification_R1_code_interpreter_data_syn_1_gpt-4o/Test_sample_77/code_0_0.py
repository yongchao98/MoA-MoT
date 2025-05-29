def apply_operations(s):
    while True:
        original = s
        # Rule 1: If the string prefix is 'ab', replace it with 'ca'.
        if s.startswith('ab'):
            s = 'ca' + s[2:]
        
        # Rule 2: If the string prefix is 'cb', replace it with 'aa' and delete the last character.
        elif s.startswith('cb'):
            s = 'aa' + s[2:-1]
        
        # Rule 3: If the string suffix is 'ac', replace it with 'cb'.
        elif s.endswith('ac'):
            s = s[:-2] + 'cb'
        
        # Rule 4: If the string prefix is 'bc', delete the first two characters and append 'aa' to the end.
        elif s.startswith('bc'):
            s = s[2:] + 'aa'
        
        # Rule 5: If the string ends with 'aa', replace it with 'cc'.
        elif s.endswith('aa'):
            s = s[:-2] + 'cc'
        
        # Rule 6: If the string starts with 'bb', remove the second character.
        elif s.startswith('bb'):
            s = s[0] + s[2:]
        
        # If no changes were made, break the loop
        if s == original:
            break
    
    return s

# Initial string
initial_string = "bbcaaccbbccbbbaa"
final_string = apply_operations(initial_string)
print(final_string)