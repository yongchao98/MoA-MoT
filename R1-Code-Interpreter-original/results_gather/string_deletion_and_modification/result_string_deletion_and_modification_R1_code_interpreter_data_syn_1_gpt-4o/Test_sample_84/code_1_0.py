def transform_string(s):
    while True:
        if 'abc' in s:
            s = s.replace('abc', 'cab', 1)
        elif s.startswith('ac'):
            s = 'zz' + s[2:]
        elif 'acb' in s:
            s = s.replace('acb', 'bca', 1)
        elif s.count('b') % 2 == 0 and s.count('b') > 0:
            s += 'ab'
        elif s.endswith('cc'):
            s = 'a' + s[:-2] + 'b'
        elif s.endswith('ba'):
            s = s[:-2] + 'ab'
        else:
            break
    return s

final_string = transform_string("abbbaaabcacbcbbcaba")
print(final_string)