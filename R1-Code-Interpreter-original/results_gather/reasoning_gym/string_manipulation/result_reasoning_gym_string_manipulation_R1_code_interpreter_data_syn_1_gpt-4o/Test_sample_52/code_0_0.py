def transform_string(s):
    seen_states = set()
    while s not in seen_states:
        seen_states.add(s)
        if 'abc' in s:
            s = s.replace('abc', 'cab', 1)
        elif 'ca' in s[1:]:
            index = s.find('ca', 1)
            s = s[:index] + s[index+2:]
        elif s.startswith('bb'):
            s = s[0] + s[2:]
        elif 'acb' in s:
            s = s.replace('acb', 'bca', 1)
        elif s.endswith('ba'):
            s = s[:-2] + 'ab'
        elif 'bca' in s:
            index = s.find('bca')
            s = s[:index] + s[index+3:]
        elif s.endswith('bb'):
            s = s[:-2]
        else:
            break
    return s

initial_string = "abaccbabbccbca"
final_string = transform_string(initial_string)
print(final_string)