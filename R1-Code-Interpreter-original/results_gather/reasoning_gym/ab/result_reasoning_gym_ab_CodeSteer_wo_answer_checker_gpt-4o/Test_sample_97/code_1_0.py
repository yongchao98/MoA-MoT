def compute_program(tokens):
    tokens = tokens.split()
    changed = True

    while changed:
        changed = False
        i = 0
        while i < len(tokens) - 1:
            if tokens[i] == 'A#' and tokens[i + 1] == '#A':
                tokens.pop(i)
                tokens.pop(i)
                changed = True
            elif tokens[i] == 'A#' and tokens[i + 1] == '#B':
                tokens[i], tokens[i + 1] = '#B', 'A#'
                i += 1
            elif tokens[i] == 'B#' and tokens[i + 1] == '#A':
                tokens[i], tokens[i + 1] = '#A', 'B#'
                i += 1
            elif tokens[i] == 'B#' and tokens[i + 1] == '#B':
                tokens.pop(i)
                tokens.pop(i)
                changed = True
            else:
                i += 1

    return ' '.join(tokens)

program = "#A A# #B B# #A B# A# B# B# #A"
final_state = compute_program(program)
print(f"<<<{final_state}>>>")