def transform_string(s):
    while True:
        initial = s
        if s.startswith('bc'):
            s = s[2:] + 'aa'
        elif s.endswith('aa'):
            s = s[:-2] + 'cc'
        elif s.startswith('aa'):
            s = s[1:]
        elif 'ca' in s[1:]:
            index = s.find('ca', 1)
            s = s[:index] + s[index+2:]
        elif s.startswith('cb'):
            s = 'aa' + s[2:-1]
        elif 'acb' in s:
            index = s.find('acb')
            s = s[:index] + 'bca' + s[index+3:]
        
        if s == initial:
            break
    return s

result = transform_string("cbcbaabacc")
print(result)