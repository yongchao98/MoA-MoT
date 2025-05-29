def transform_string(s):
    while True:
        if 'bca' in s:
            s = s.replace('bca', '', 1)
        elif s.endswith('bb'):
            s = s[:-2]
        elif s.endswith('ba'):
            s = s[:-2] + 'ab'
        elif s.endswith('ca'):
            s = s[:-1]
        elif 'acb' in s:
            s = s.replace('acb', 'bca', 1)
        elif s.startswith('ac'):
            s = 'zz' + s[2:]
        else:
            break
    return s

initial_string = "acacabababcacbcc"
final_string = transform_string(initial_string)
print(f'<<< "{final_string}" >>>')