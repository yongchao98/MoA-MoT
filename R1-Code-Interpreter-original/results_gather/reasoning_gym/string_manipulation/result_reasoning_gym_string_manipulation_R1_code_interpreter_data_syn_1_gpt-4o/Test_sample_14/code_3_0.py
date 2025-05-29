def transform_string(s):
    seen_states = set()
    while s not in seen_states:
        seen_states.add(s)
        original_s = s
        
        # Rule 1
        if len(s) > 15:
            mid_index = len(s) // 2
            s = s[:mid_index] + s[mid_index+1:]
        
        # Rule 2
        elif s.count('b') % 2 == 0 and s.count('b') > 0:
            s += 'ab'
        
        # Rule 3
        elif s.endswith('ac'):
            s = s[:-2] + 'cb'
        
        # Rule 4
        elif s.startswith('ab'):
            s = 'ca' + s[2:]
        
        # If no rule was applied, break
        if s == original_s:
            break
    
    return s

final_string = transform_string("baccacc")
print(final_string)