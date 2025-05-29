def process_string(s):
    while True:
        if s.startswith('aa'):
            s = s[1:]
        elif s.count('b') % 2 == 0 and 'b' in s:
            s += 'ab'
        elif s.startswith('ab'):
            s = 'ca' + s[2:]
        elif 'acb' in s:
            s = s.replace('acb', 'bca', 1)
        elif 'ca' in s[1:]:
            index = s.find('ca', 1)
            s = s[:index] + s[index+2:]
        elif s.startswith('cb'):
            s = 'aa' + s[2:-1]
        else:
            break
    return s

initial_string = "acbbbcbbaaabbbc"
final_string = process_string(initial_string)
print(f'<<< "{final_string}" >>>')