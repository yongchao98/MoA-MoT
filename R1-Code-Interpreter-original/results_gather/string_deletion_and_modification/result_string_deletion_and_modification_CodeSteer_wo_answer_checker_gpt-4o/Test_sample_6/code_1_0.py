def transform_string(s):
    while True:
        original = s
        if s.startswith('bc'):
            s = s[2:] + 'aa'
        elif 'abc' in s:
            s = s.replace('abc', 'cab', 1)
        elif s.endswith('ac'):
            s = s[:-2] + 'cb'
        elif s.startswith('ab'):
            s = 'ca' + s[2:]
        elif 'acb' in s:
            s = s.replace('acb', 'bca', 1)
        elif s.endswith('bb'):
            s = s[:-2]
        
        # If no changes were made, break the loop
        if s == original:
            break
    
    return s

initial_string = "acbbaacabcacbbccaaa"
final_string = transform_string(initial_string)
print(f'<<< "{final_string}" >>>')