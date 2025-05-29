def compute_program(program):
    tokens = program.split()
    changed = True
    while changed:
        changed = False
        i = 0
        while i < len(tokens) - 1:
            if tokens[i] == 'A#' and tokens[i + 1] == '#A':
                # Apply rule: A# #A -> nothing
                tokens.pop(i)
                tokens.pop(i)
                changed = True
            elif tokens[i] == 'A#' and tokens[i + 1] == '#B':
                # Apply rule: A# #B -> #B A#
                tokens[i], tokens[i + 1] = '#B', 'A#'
                changed = True
                i += 1
            elif tokens[i] == 'B#' and tokens[i + 1] == '#A':
                # Apply rule: B# #A -> #A B#
                tokens[i], tokens[i + 1] = '#A', 'B#'
                changed = True
                i += 1
            elif tokens[i] == 'B#' and tokens[i + 1] == '#B':
                # Apply rule: B# #B -> nothing
                tokens.pop(i)
                tokens.pop(i)
                changed = True
            else:
                i += 1
    return ' '.join(tokens)

program = "#B #A #B #B A# #A B# #A B# #A"
final_state = compute_program(program)
print(final_state)