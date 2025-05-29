def transform_string(s):
    previous_states = set()
    while s not in previous_states:
        previous_states.add(s)
        if s.endswith('ac'):
            s = s[:-2] + 'cb'
        elif s.endswith('bb'):
            s = s[:-2]
        elif s.endswith('ba'):
            s = s[:-2] + 'ab'
        elif 'ca' in s[1:]:
            index = s.find('ca', 1)
            s = s[:index] + s[index+2:]
        elif s.count('b') % 2 == 0 and s.count('b') > 0:
            s += 'ab'
        elif s.endswith('ca'):
            s = s[:-1]
        elif s.startswith('ab'):
            s = 'ca' + s[2:]
        elif 'abc' in s:
            index = s.find('abc')
            s = s[:index] + 'cab' + s[index+3:]
        else:
            break
    return s

result = transform_string("aabcaba")
print(result)