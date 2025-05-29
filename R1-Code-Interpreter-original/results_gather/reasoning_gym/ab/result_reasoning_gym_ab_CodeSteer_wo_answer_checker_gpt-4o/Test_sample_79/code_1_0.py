def compute_program(program):
    tokens = program.split()
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
                tokens[i], tokens[i + 1] = '#B', 'A#'
                i += 1
                changed = True
            elif tokens[i] == 'B#' and tokens[i + 1] == '#A':
                # B# #A becomes #A B#
                tokens[i], tokens[i + 1] = '#A', 'B#'
                i += 1
                changed = True
            elif tokens[i] == 'B#' and tokens[i + 1] == '#B':
                # B# #B becomes nothing
                tokens.pop(i)
                tokens.pop(i)
                changed = True
            else:
                i += 1

    return ' '.join(tokens)

program = "B# #A #B B# #B #B B# B# #B B#"
final_state = compute_program(program)
print(f"<<<{final_state}>>>")