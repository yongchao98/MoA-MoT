def transform_string(s):
    previous_states = set()
    
    while s not in previous_states:
        previous_states.add(s)
        
        if s.startswith('bc'):
            s = s[2:] + 'aa'
        elif s.endswith('ba'):
            s = s[:-2] + 'ab'
        elif s.endswith('ca'):
            s = s[:-1]
        elif s.startswith('ac'):
            s = 'zz' + s[2:]
        else:
            break
    
    return s

initial_string = "bbcacbbacbaccaababba"
final_string = transform_string(initial_string)
print(final_string)