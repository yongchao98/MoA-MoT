def transform_string(s):
    previous_states = set()
    while s not in previous_states:
        previous_states.add(s)
        if s.startswith('bb'):
            s = s[0] + s[2:]  # Remove the second character
        elif 'ca' in s[1:]:
            index = s.find('ca', 1)
            s = s[:index] + s[index+2:]  # Remove 'ca'
        elif s.endswith('bb'):
            s = s[:-2]  # Remove the last two characters
        elif s.startswith('ac'):
            s = 'zz' + s[2:]  # Replace 'ac' with 'zz'
        elif 'bca' in s:
            index = s.find('bca')
            s = s[:index] + s[index+3:]  # Remove 'bca'
        else:
            break
    return s

final_string = transform_string("cccccac")
print(final_string)