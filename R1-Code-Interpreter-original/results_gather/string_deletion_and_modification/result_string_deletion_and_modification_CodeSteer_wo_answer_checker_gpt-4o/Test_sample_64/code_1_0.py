def process_string(s):
    while True:
        if s.endswith('cc'):
            s = 'a' + s[:-2] + 'b'
        elif len(s) > 15:
            mid_index = len(s) // 2
            s = s[:mid_index] + s[mid_index+1:]
        elif 'acb' in s:
            s = s.replace('acb', 'bca', 1)
        elif s.startswith('cb'):
            s = 'aa' + s[2:-1]
        elif s.startswith('bb'):
            s = s[0] + s[2:]
        elif s.startswith('aa'):
            s = s[1:]
        else:
            break
    return s

initial_string = "baccabcaabcaabbcbca"
final_string = process_string(initial_string)
print(f'<<< "{final_string}" >>>')