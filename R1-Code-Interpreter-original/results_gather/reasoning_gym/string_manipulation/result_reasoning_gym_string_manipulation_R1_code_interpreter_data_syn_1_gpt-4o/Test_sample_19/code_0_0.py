def transform_string(s):
    seen_states = set()
    while True:
        if s in seen_states:
            break
        seen_states.add(s)
        
        original = s
        
        # Rule 1
        if 'ca' in s[1:]:
            index = s[1:].index('ca') + 1
            s = s[:index] + s[index+2:]
            continue
        
        # Rule 2
        if s.endswith('ba'):
            s = s[:-2] + 'ab'
            continue
        
        # Rule 3
        if s.startswith('cc'):
            s = s[2:]
            continue
        
        # Rule 4
        if 'acb' in s:
            index = s.index('acb')
            s = s[:index] + 'bca' + s[index+3:]
            continue
        
        # Rule 5
        if 'abc' in s:
            index = s.index('abc')
            s = s[:index] + 'cab' + s[index+3:]
            continue
        
        # Rule 6
        if s.endswith('cc'):
            s = 'a' + s[:-2] + 'b'
            continue
        
        # Rule 7
        if s.count('b') % 2 == 0 and s.count('b') > 0:
            s += 'ab'
            continue
        
        # If no rule applied, break
        if s == original:
            break
    
    return s

final_string = transform_string("ccaaacbbcbcc")
print(final_string)