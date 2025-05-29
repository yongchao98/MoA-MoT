def transform_string(s):
    previous_states = set()
    while s not in previous_states:
        previous_states.add(s)
        if 'bca' in s:
            s = s.replace('bca', '', 1)
        elif s.endswith('ca'):
            s = s[:-1]
        elif s.endswith('ac'):
            s = s[:-2] + 'cb'
        elif 'ca' in s[1:]:
            index = s.find('ca', 1)
            s = s[:index] + s[index+2:]
        elif s.startswith('ca'):
            s = 'bb' + s[2:] + 'c'
        elif s.startswith('cc'):
            s = s[2:]
        else:
            break
    return s

final_string = transform_string("ccbacccca")
print(final_string)