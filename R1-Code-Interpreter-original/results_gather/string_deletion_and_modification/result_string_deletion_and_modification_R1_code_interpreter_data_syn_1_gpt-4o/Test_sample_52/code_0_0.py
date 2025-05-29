def transform_string(s):
    while True:
        if s.endswith('bb'):
            s = s[:-2]
        elif s.endswith('ac'):
            s = s[:-2] + 'cb'
        elif 'acb' in s:
            s = s.replace('acb', 'bca', 1)
        elif s.startswith('bc'):
            s = s[2:] + 'aa'
        elif s.startswith('ab'):
            s = 'ca' + s[2:]
        elif s.startswith('ca'):
            s = 'bb' + s[2:] + 'c'
        else:
            break
    return s

initial_string = "bccaaaacaa"
final_string = transform_string(initial_string)
print(final_string)