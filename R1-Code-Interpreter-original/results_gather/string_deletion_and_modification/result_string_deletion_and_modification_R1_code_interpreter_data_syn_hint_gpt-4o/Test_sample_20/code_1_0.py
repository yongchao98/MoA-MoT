def process_string(s):
    while True:
        if 'bca' in s:
            s = s.replace('bca', '', 1)
        elif s.count('b') % 2 == 0 and s.count('b') > 0:
            s += 'ab'
        elif s.startswith('aa'):
            s = s[1:]
        elif s.endswith('aa'):
            s = s[:-2] + 'cc'
        elif s.startswith('bb'):
            s = s[0] + s[2:]
        elif 'acb' in s:
            s = s.replace('acb', 'bca', 1)
        else:
            break
    return s

final_string = process_string("baacbabbabbbc")
print(final_string)