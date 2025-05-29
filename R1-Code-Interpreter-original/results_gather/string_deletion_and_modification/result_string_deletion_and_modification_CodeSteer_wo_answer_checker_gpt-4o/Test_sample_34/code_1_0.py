def transform_string(s):
    while True:
        original = s
        # Operation 1
        if s.startswith('cb'):
            s = 'aa' + s[2:-1]
        # Operation 2
        elif s.endswith('bb'):
            s = s[:-2]
        # Operation 3
        elif 'acb' in s:
            s = s.replace('acb', 'bca', 1)
        # Operation 4
        elif s.startswith('bc'):
            s = s[2:] + 'aa'
        # Operation 5
        elif s.startswith('ac'):
            s = 'zz' + s[2:]
        # Operation 6
        elif 'ca' in s[1:]:
            index = s.find('ca', 1)
            s = s[:index] + s[index+2:]
        # If no operations were applied, break the loop
        if s == original:
            break
    return s

final_string = transform_string("aacccabbbbcaacaacbcc")
print(f'<<< "{final_string}" >>>')