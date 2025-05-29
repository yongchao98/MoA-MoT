def transform_string(s):
    while True:
        if s.endswith('bb'):
            s = s[:-2]
        elif 'ca' in s[1:]:
            s = s[:s.index('ca', 1)] + s[s.index('ca', 1) + 2:]
        elif s.endswith('cc'):
            s = 'a' + s[:-2] + 'b'
        elif s.endswith('ac'):
            s = s[:-2] + 'cb'
        elif s.endswith('ca'):
            s = s[:-1]
        elif s.startswith('cb'):
            s = 'aa' + s[2:-1]
        else:
            break
    return s

initial_string = "bbaaacccaaccabbca"
final_string = transform_string(initial_string)
print(f'<<< "{final_string}" >>>')