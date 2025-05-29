def transform_string(s):
    while True:
        original = s
        if s.startswith('ac'):
            s = 'zz' + s[2:]
        elif 'ca' in s[1:]:
            index = s.find('ca', 1)
            s = s[:index] + s[index+2:]
        elif 'bca' in s:
            index = s.find('bca')
            s = s[:index] + s[index+3:]
        elif s.endswith('bb'):
            s = s[:-2]
        elif s.endswith('aa'):
            s = s[:-2] + 'cc'
        elif s.startswith('cc'):
            s = s[2:]
        
        if s == original:
            break
    return s

final_string = transform_string("acbbabbcabccacaccc")
print(f'<<< "{final_string}" >>>')