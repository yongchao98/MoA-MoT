def transform_string(s):
    previous_states = set()
    
    while s not in previous_states:
        previous_states.add(s)
        
        if s.endswith('ba'):
            s = s[:-2] + 'ab'
        elif 'ca' in s[1:]:
            index = s.find('ca', 1)
            s = s[:index] + s[index+2:]
        elif s.endswith('aa'):
            s = s[:-2] + 'cc'
        elif s.startswith('bb'):
            s = s[0] + s[2:]
        elif 'bca' in s:
            index = s.find('bca')
            s = s[:index] + s[index+3:]
        elif s.startswith('ac'):
            s = 'zz' + s[2:]
        elif len(s) > 15:
            mid = len(s) // 2
            s = s[:mid] + s[mid+1:]
        else:
            break
    
    return s

final_string = transform_string("cabbcccacc")
print(final_string)