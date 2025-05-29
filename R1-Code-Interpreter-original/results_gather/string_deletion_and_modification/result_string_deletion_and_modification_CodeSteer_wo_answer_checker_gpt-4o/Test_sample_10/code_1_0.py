def transform_string(s):
    while True:
        if s.endswith('ba'):
            s = s[:-2] + 'ab'
        elif s.startswith('aa'):
            s = s[1:]
        elif s.endswith('ac'):
            s = s[:-2] + 'cb'
        elif s.startswith('cc'):
            s = s[2:]
        elif s.startswith('bc'):
            s = s[2:] + 'aa'
        elif s.endswith('cc'):
            s = 'a' + s[:-2] + 'b'
        else:
            break
    return s

initial_string = "ccccaaaaccbabcbabcbb"
final_string = transform_string(initial_string)
print(f'<<< "{final_string}" >>>')