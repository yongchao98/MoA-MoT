def transform_string(s):
    while True:
        if s.endswith('ac'):
            s = s[:-2] + 'cb'
        elif s.startswith('bc'):
            s = s[2:] + 'aa'
        elif 'abc' in s:
            s = s.replace('abc', 'cab', 1)
        elif s.endswith('cc'):
            s = 'a' + s[:-2] + 'b'
        elif s.endswith('ca'):
            s = s[:-1]
        elif 'ca' in s[1:]:
            index = s.find('ca', 1)
            s = s[:index] + s[index+2:]
        else:
            break
    return s

initial_string = "accacabbccbbbacaabc"
final_string = transform_string(initial_string)
print(final_string)