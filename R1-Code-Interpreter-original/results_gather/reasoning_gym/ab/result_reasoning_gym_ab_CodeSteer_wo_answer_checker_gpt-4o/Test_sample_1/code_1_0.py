def compute_program(program):
    tokens = program.split()
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
                changed = True
                i += 1
            elif tokens[i] == 'B#' and tokens[i + 1] == '#A':
                tokens[i], tokens[i + 1] = '#A', 'B#'
                changed = True
                i += 1
            elif tokens[i] == 'B#' and tokens[i + 1] == '#B':
                tokens.pop(i)
                tokens.pop(i)
                changed = True
            else:
                i += 1

    return ' '.join(tokens)

program = "B# B# #A #B B# #A #B #B #A A#"
final_state = compute_program(program)
print(f"<<<{final_state}>>>")