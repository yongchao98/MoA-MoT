def compute_program(program):
    tokens = program.split()
    changed = True

    while changed:
        changed = False
        i = 0
        while i < len(tokens) - 1:
            pair = (tokens[i], tokens[i + 1])
            if pair == ('A#', '#A'):
                tokens.pop(i)
                tokens.pop(i)
                changed = True
            elif pair == ('A#', '#B'):
                tokens[i], tokens[i + 1] = '#B', 'A#'
                i += 1
                changed = True
            elif pair == ('B#', '#A'):
                tokens[i], tokens[i + 1] = '#A', 'B#'
                i += 1
                changed = True
            elif pair == ('B#', '#B'):
                tokens.pop(i)
                tokens.pop(i)
                changed = True
            else:
                i += 1

    return ' '.join(tokens)

program = "A# A# #A B# B# B# A# A# #B A#"
final_state = compute_program(program)
print(f"<<<{final_state}>>>")