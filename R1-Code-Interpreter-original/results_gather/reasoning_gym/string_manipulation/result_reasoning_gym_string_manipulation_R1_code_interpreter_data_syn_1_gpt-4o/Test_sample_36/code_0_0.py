def transform_string(s):
    previous_states = set()
    while True:
        if s in previous_states:
            break
        previous_states.add(s)
        
        # Apply Rule 1
        if s.endswith('ac'):
            s = s[:-2] + 'cb'
            continue
        
        # Apply Rule 2
        if s.count('b') % 2 == 0 and s.count('b') > 0:
            s += 'ab'
            continue
        
        # Apply Rule 3
        if s.endswith('aa'):
            s = s[:-2] + 'cc'
            continue
        
        # If no rule was applied, break the loop
        break
    
    return s

initial_string = "cccaccabbbccbbcbacc"
final_string = transform_string(initial_string)
print(final_string)