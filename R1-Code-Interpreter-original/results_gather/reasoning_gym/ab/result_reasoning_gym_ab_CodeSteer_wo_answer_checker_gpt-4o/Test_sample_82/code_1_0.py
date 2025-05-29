def compute_program(tokens):
    tokens = tokens.split()
    changed = True

    while changed:
        changed = False
        i = 0
        while i < len(tokens) - 1:
            if tokens[i] == 'A#' and tokens[i + 1] == '#A':
                # A# #A becomes nothing
                tokens.pop(i)
                tokens.pop(i)
                changed = True
            elif tokens[i] == 'A#' and tokens[i + 1] == '#B':
                # A# #B becomes #B A#
                tokens[i], tokens[i + 1] = tokens[i + 1], tokens[i]
                i += 1
            elif tokens[i] == 'B#' and tokens[i + 1] == '#A':
                # B# #A becomes #A B#
                tokens[i], tokens[i + 1] = tokens[i + 1], tokens[i]
                i += 1
            elif tokens[i] == 'B#' and tokens[i + 1] == '#B':
                # B# #B becomes nothing
                tokens.pop(i)
                tokens.pop(i)
                changed = True
            else:
                i += 1

    return ' '.join(tokens)

program = "B# A# #A B# A# B# #A A# #A #A"
final_state = compute_program(program)
print(f"<<<{final_state}>>>")