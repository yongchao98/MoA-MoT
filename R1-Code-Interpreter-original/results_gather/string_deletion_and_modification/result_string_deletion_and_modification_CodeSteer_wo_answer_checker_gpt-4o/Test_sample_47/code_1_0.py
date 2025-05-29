def transform_string(s):
    while True:
        original = s
        if s.startswith('aa'):
            s = s[1:]
        elif s.startswith('cb'):
            s = 'aa' + s[2:-1]
        elif s.startswith('cc'):
            s = s[2:]
        elif s.endswith('bb'):
            s = s[:-2]
        elif 'acb' in s:
            s = s.replace('acb', 'bca', 1)
        elif s.endswith('aa'):
            s = s[:-2] + 'cc'
        
        # If no changes were made, break the loop
        if s == original:
            break
    return s

final_string = transform_string("cacbcbbbabbabacc")
print(f'<<< "{final_string}" >>>')