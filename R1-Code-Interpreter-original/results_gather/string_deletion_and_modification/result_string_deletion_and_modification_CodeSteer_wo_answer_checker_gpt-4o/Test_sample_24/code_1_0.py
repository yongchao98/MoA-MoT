def transform_string(s):
    while True:
        if s.endswith('ba'):
            s = s[:-2] + 'ab'
        elif 'acb' in s:
            s = s.replace('acb', 'bca', 1)
        elif 'bca' in s:
            s = s.replace('bca', '', 1)
        elif s.startswith('ab'):
            s = 'ca' + s[2:]
        elif s.startswith('bb'):
            s = s[0] + s[2:]
        elif s.startswith('ac'):
            s = 'zz' + s[2:]
        else:
            break
    return s

final_string = transform_string("acbacbcaacabcababbc")
print(f'<<< "{final_string}" >>>')