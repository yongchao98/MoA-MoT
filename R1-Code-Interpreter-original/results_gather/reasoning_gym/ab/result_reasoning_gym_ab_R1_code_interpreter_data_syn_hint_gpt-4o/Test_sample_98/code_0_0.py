def compute_program(program):
    tokens = program.split()
    i = 0
    while i < len(tokens) - 1:
        if tokens[i] == 'A#' and tokens[i + 1] == '#A':
            # Apply rule: A# #A -> nothing
            tokens.pop(i)
            tokens.pop(i)
            i = max(i - 1, 0)  # Move back one step to recheck
        elif tokens[i] == 'A#' and tokens[i + 1] == '#B':
            # Apply rule: A# #B -> #B A#
            tokens[i], tokens[i + 1] = '#B', 'A#'
            i += 1
        elif tokens[i] == 'B#' and tokens[i + 1] == '#A':
            # Apply rule: B# #A -> #A B#
            tokens[i], tokens[i + 1] = '#A', 'B#'
            i += 1
        elif tokens[i] == 'B#' and tokens[i + 1] == '#B':
            # Apply rule: B# #B -> nothing
            tokens.pop(i)
            tokens.pop(i)
            i = max(i - 1, 0)  # Move back one step to recheck
        else:
            i += 1
    return ' '.join(tokens)

program = "#B #A #B #B A# #A B# #A B# #A"
final_state = compute_program(program)
print(final_state)