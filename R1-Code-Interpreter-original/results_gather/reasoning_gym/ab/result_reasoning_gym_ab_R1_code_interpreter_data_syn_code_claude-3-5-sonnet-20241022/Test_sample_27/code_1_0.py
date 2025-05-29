def apply_rules(program):
    tokens = program.split()
    i = 0
    while i < len(tokens) - 1:
        pair = tokens[i] + tokens[i+1]
        if 'A##A' in pair:  # A# #A
            tokens.pop(i)
            tokens.pop(i)
            i = max(0, i-1)  # Go back one step to check new adjacencies
        elif 'A##B' in pair:  # A# #B
            tokens[i:i+2] = ['#B', 'A#']
            i = max(0, i-1)
        elif 'B##A' in pair:  # B# #A
            tokens[i:i+2] = ['#A', 'B#']
            i = max(0, i-1)
        elif 'B##B' in pair:  # B# #B
            tokens.pop(i)
            tokens.pop(i)
            i = max(0, i-1)
        else:
            i += 1
    return ' '.join(tokens)

program = "B# A# #B A# B# #A A# A# #B #B"
print(apply_rules(program))