def transform_string(s):
    while True:
        original = s
        if 'acb' in s:
            s = s.replace('acb', 'bca', 1)
        elif s.startswith('aa'):
            s = s[1:]
        elif s.startswith('ca'):
            s = 'bb' + s[2:] + 'c'
        elif s.startswith('bb'):
            s = s[0] + s[2:]
        elif s.endswith('bb'):
            s = s[:-2]
        elif s.startswith('bc'):
            s = s[2:] + 'aa'
        
        # If no changes were made, break the loop
        if s == original:
            break
    
    return s

initial_string = "ccacbaabacbcacbccc"
final_string = transform_string(initial_string)
print(f'<<< "{final_string}" >>>')