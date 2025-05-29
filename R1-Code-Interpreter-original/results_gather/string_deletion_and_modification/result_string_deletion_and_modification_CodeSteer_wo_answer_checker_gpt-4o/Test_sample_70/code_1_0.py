def transform_string(s):
    while True:
        if 'ca' in s[1:]:
            s = s[:s.index('ca', 1)] + s[s.index('ca', 1) + 2:]
        elif 'acb' in s:
            s = s.replace('acb', 'bca', 1)
        elif 'abc' in s:
            s = s.replace('abc', 'cab', 1)
        elif s.startswith('bc'):
            s = s[2:] + 'aa'
        elif s.startswith('aa'):
            s = s[1:]
        elif s.startswith('cc'):
            s = s[2:]
        else:
            break
    return s

initial_string = "acaabacaccbcbccbcaa"
final_string = transform_string(initial_string)
print(f'<<< "{final_string}" >>>')