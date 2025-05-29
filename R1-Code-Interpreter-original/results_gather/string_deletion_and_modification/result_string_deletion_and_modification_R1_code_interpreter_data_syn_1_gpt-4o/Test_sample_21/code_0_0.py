def transform_string(s):
    while True:
        if s.endswith('aa'):
            s = s[:-2] + 'cc'
        elif s.startswith('cc'):
            s = s[2:]
        elif 'abc' in s:
            s = s.replace('abc', 'cab', 1)
        elif 'bca' in s:
            s = s.replace('bca', '', 1)
        elif s.startswith('ac'):
            s = 'zz' + s[2:]
        elif 'ca' in s[1:]:
            index = s[1:].index('ca') + 1
            s = s[:index] + s[index+2:]
        else:
            break
    return s

initial_string = "acacbccaabaacbbaa"
final_string = transform_string(initial_string)
print(final_string)