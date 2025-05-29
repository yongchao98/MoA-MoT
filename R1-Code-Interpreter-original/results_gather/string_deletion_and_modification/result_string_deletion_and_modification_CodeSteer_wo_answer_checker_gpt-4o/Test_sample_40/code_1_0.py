def transform_string(s):
    while True:
        original = s
        if 'abc' in s:
            s = s.replace('abc', 'cab', 1)
        elif 'ca' in s[1:]:
            index = s.find('ca', 1)
            s = s[:index] + s[index+2:]
        elif s.endswith('ca'):
            s = s[:-1]
        elif s.startswith('bb'):
            s = s[0] + s[2:]
        elif s.startswith('ca'):
            s = 'bb' + s[2:] + 'c'
        elif s.endswith('ac'):
            s = s[:-2] + 'cb'
        
        if s == original:
            break
    return s

final_string = transform_string("bbaacaabcbbba")
print(f'<<< "{final_string}" >>>')