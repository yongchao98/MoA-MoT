def compute_program(program):
    tokens = program.split()
    i = 0
    while i < len(tokens) - 1:
        if tokens[i] == 'A#' and tokens[i + 1] == '#A':
            tokens.pop(i)
            tokens.pop(i)
            i = 0  # Restart from the beginning
        elif tokens[i] == 'A#' and tokens[i + 1] == '#B':
            tokens[i], tokens[i + 1] = '#B', 'A#'
            i = 0  # Restart from the beginning
        elif tokens[i] == 'B#' and tokens[i + 1] == '#A':
            tokens[i], tokens[i + 1] = '#A', 'B#'
            i = 0  # Restart from the beginning
        elif tokens[i] == 'B#' and tokens[i + 1] == '#B':
            tokens.pop(i)
            tokens.pop(i)
            i = 0  # Restart from the beginning
        else:
            i += 1
    return ' '.join(tokens)

program = "B# #B #B #A B# #B A# #A #A A#"
final_state = compute_program(program)
print(final_state)