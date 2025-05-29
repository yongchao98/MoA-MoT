def compute_program(program):
    tokens = program.split()
    i = 0
    while i < len(tokens) - 1:
        if tokens[i] == 'A#' and tokens[i + 1] == '#A':
            # Apply rule: A# #A becomes nothing
            tokens.pop(i)
            tokens.pop(i)
            i = max(i - 1, 0)  # Step back to check for new possible rewrites
        elif tokens[i] == 'A#' and tokens[i + 1] == '#B':
            # Apply rule: A# #B becomes #B A#
            tokens[i], tokens[i + 1] = '#B', 'A#'
            i += 1
        elif tokens[i] == 'B#' and tokens[i + 1] == '#A':
            # Apply rule: B# #A becomes #A B#
            tokens[i], tokens[i + 1] = '#A', 'B#'
            i += 1
        elif tokens[i] == 'B#' and tokens[i + 1] == '#B':
            # Apply rule: B# #B becomes nothing
            tokens.pop(i)
            tokens.pop(i)
            i = max(i - 1, 0)  # Step back to check for new possible rewrites
        else:
            i += 1
    return ' '.join(tokens)

program = "A# A# B# A# #A #A #B #B #B #B"
final_state = compute_program(program)
print(final_state)