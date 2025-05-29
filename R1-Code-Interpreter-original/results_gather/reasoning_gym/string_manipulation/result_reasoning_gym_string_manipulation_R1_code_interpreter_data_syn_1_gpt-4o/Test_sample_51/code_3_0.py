def transform_string(s):
    seen_states = set()
    while s not in seen_states:
        seen_states.add(s)
        print(f"Current string: {s}")
        if s.startswith('aa'):
            print(f"Rule 1 applied on {s}")
            s = s[1:]
        elif s.startswith('bb'):
            print(f"Rule 2 applied on {s}")
            s = s[0] + s[2:]
        elif s.startswith('ca'):
            print(f"Rule 3 applied on {s}")
            s = 'bb' + s[2:] + 'c'
        elif s.endswith('ac'):
            print(f"Rule 4 applied on {s}")
            s = s[:-2] + 'cb'
        elif s.startswith('ab'):
            print(f"Rule 5 applied on {s}")
            s = 'ca' + s[2:]
        elif s.endswith('aa'):
            print(f"Rule 6 applied on {s}")
            s = s[:-2] + 'cc'
        elif 'acb' in s:
            print(f"Rule 7 applied on {s}")
            s = s.replace('acb', 'bca', 1)
        elif s.startswith('ac'):
            print(f"Rule 8 applied on {s}")
            s = 'zz' + s[2:]
        else:
            print("No applicable rule found.")
            break
    return s

initial_string = "baabaacaaabaababbcab"
final_string = transform_string(initial_string)
print(final_string)