def transform_string(s):
    while True:
        if s.startswith('ab'):
            s = 'ca' + s[2:]
        elif len(s) > 15:
            mid_index = len(s) // 2
            s = s[:mid_index] + s[mid_index+1:]
        elif s.endswith('aa'):
            s = s[:-2] + 'cc'
        elif 'ca' in s[1:]:
            first_ca_index = s[1:].index('ca') + 1
            s = s[:first_ca_index] + s[first_ca_index+2:]
        elif s.count('b') % 2 == 0 and s.count('b') > 0:
            s += 'ab'
        elif s.endswith('cc'):
            s = 'a' + s[:-2] + 'b'
        else:
            break
    return s

final_string = transform_string("cabacbbababcbacb")
print(f'<<< "{final_string}" >>>')