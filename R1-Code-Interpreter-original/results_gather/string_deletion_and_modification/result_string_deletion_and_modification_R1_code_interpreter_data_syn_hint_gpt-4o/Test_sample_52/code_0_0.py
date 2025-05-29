def transform_string(s):
    while True:
        if s.startswith('bc'):
            s = s[2:] + 'aa'
        elif s.startswith('ca'):
            s = 'bb' + s[2:] + 'c'
        elif s.startswith('ab'):
            s = 'ca' + s[2:]
        elif s.endswith('bb'):
            s = s[:-2]
        elif s.endswith('ac'):
            s = s[:-2] + 'cb'
        elif 'acb' in s:
            s = s.replace('acb', 'bca', 1)
        else:
            break
    return s

final_string = transform_string("bccaaaacaa")
print(final_string)