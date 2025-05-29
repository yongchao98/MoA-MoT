def compute_program(tokens):
    i = 0
    while i < len(tokens) - 1:
        if tokens[i] == 'A#' and tokens[i + 1] == '#A':
            # Apply rule: A# #A -> nothing
            tokens.pop(i)
            tokens.pop(i)
            i = max(i - 1, 0)  # Move back one step to check for new pairs
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
            i = max(i - 1, 0)  # Move back one step to check for new pairs
        else:
            i += 1
    return tokens

# Initial program
program = ['#A', 'A#', '#B', '#A', '#B', 'A#', 'B#', 'B#', '#B', 'B#']
final_state = compute_program(program)
print(final_state)