def transform_string(s):
    while True:
        if s.endswith('ac'):
            s = s[:-2] + 'cb'
        elif s.startswith('bc'):
            s = s[2:] + 'aa'
        elif s.endswith('bb'):
            s = s[:-2]
        elif s.startswith('ab'):
            s = 'ca' + s[2:]
        elif s.endswith('cc'):
            s = 'a' + s[:-2] + 'b'
        elif s.startswith('cb'):
            s = 'aa' + s[2:-1]
        else:
            break
    return s

initial_string = "cbcabbabcccaaaacbb"
final_string = transform_string(initial_string)
print(f'<<< "{final_string}" >>>')