def transform_string(s):
    previous_states = set()
    
    while s not in previous_states:
        previous_states.add(s)
        
        if s.endswith('ac'):
            s = s[:-2] + 'cb'
        elif 'abc' in s:
            s = s.replace('abc', 'cab', 1)
        elif s.startswith('cb'):
            s = 'aa' + s[2:-1]
        elif s.startswith('bb'):
            s = s[0] + s[2:]
        elif s.startswith('ac'):
            s = 'zz' + s[2:]
        elif s.endswith('ba'):
            s = s[:-2] + 'ab'
        elif len(s) > 15:
            mid = len(s) // 2
            s = s[:mid] + s[mid+1:]
        elif s.endswith('ca'):
            s = s[:-1]
        else:
            break
    
    return s

final_string = transform_string("bbabcca")
print(final_string)