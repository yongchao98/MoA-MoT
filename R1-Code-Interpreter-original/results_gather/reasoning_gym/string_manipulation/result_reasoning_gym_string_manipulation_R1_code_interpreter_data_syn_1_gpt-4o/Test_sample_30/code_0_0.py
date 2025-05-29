def transform_string(s):
    seen_states = set()
    while True:
        if s in seen_states:
            break
        seen_states.add(s)
        
        if 'abc' in s:
            s = s.replace('abc', 'cab', 1)
        elif s.startswith('ca'):
            s = 'bb' + s[2:] + 'c'
        elif s.endswith('bb'):
            s = s[:-2]
        else:
            break
    return s

initial_string = "aabbcbacaa"
final_string = transform_string(initial_string)
print(final_string)