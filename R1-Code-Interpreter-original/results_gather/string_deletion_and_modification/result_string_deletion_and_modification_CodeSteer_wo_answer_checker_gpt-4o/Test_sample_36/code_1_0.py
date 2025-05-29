def process_string(s):
    while True:
        original = s
        if 'acb' in s:
            s = s.replace('acb', 'bca', 1)
        elif 'bca' in s:
            s = s.replace('bca', '', 1)
        elif len(s) > 15:
            mid = len(s) // 2
            s = s[:mid] + s[mid+1:]
        elif s.startswith('ca'):
            s = 'bb' + s[2:] + 'c'
        elif s.endswith('bb'):
            s = s[:-2]
        elif s.startswith('bc'):
            s = s[2:] + 'aa'
        else:
            break
    return s

final_string = process_string("bcaccbaaabaacacbb")
print(f'<<< "{final_string}" >>>')