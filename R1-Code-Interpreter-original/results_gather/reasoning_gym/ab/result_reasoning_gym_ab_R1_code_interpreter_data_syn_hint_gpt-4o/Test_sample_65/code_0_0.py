def compute_program(program):
    tokens = program.split()
    i = 0
    while i < len(tokens) - 1:
        if tokens[i] == 'A#' and tokens[i + 1] == '#A':
            # Apply rule: A# #A becomes nothing
            tokens.pop(i)
            tokens.pop(i)
            i = 0  # Restart from the beginning
        elif tokens[i] == 'A#' and tokens[i + 1] == '#B':
            # Apply rule: A# #B becomes #B A#
            tokens[i], tokens[i + 1] = '#B', 'A#'
            i = 0  # Restart from the beginning
        elif tokens[i] == 'B#' and tokens[i + 1] == '#A':
            # Apply rule: B# #A becomes #A B#
            tokens[i], tokens[i + 1] = '#A', 'B#'
            i = 0  # Restart from the beginning
        elif tokens[i] == 'B#' and tokens[i + 1] == '#B':
            # Apply rule: B# #B becomes nothing
            tokens.pop(i)
            tokens.pop(i)
            i = 0  # Restart from the beginning
        else:
            i += 1  # Move to the next pair
    return ' '.join(tokens)

program = "B# #B #A A# A# A# A# #A #A A#"
final_state = compute_program(program)
print(final_state)