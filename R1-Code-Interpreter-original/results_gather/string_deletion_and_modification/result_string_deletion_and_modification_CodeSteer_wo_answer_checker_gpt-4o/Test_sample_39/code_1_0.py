def transform_string(s):
    while True:
        if s.startswith('bc'):
            s = s[2:] + 'aa'
        elif 'bca' in s:
            s = s.replace('bca', '', 1)
        elif s.startswith('ac'):
            s = 'zz' + s[2:]
        elif 'abc' in s:
            s = s.replace('abc', 'cab', 1)
        elif s.startswith('ab'):
            s = 'ca' + s[2:]
        elif s.endswith('ac'):
            s = s[:-2] + 'cb'
        else:
            break
    return s

initial_string = "bcbababbcacbbabaaac"
final_string = transform_string(initial_string)
print(f'<<< "{final_string}" >>>')