def transform_string(s):
    while True:
        if 'ca' in s[1:]:
            s = s[:s.index('ca', 1)] + s[s.index('ca', 1) + 2:]
        elif 'acb' in s:
            s = s.replace('acb', 'bca', 1)
        elif s.endswith('aa'):
            s = s[:-2] + 'cc'
        elif s.count('b') % 2 == 0 and s.count('b') > 0:
            s += 'ab'
        elif s.startswith('ca'):
            s = 'bb' + s[2:] + 'c'
        elif s.startswith('cc'):
            s = s[2:]
        else:
            break
    return s

final_string = transform_string("acbbcabcccbcc")
print(f'<<< "{final_string}" >>>')