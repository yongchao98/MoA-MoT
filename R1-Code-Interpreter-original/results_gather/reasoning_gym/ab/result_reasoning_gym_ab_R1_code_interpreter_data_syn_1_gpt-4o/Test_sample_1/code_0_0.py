def compute_program(tokens):
    i = 0
    while i < len(tokens) - 1:
        if tokens[i] == 'A#' and tokens[i + 1] == '#A':
            # A# #A becomes nothing
            tokens.pop(i)
            tokens.pop(i)
            i = max(i - 1, 0)  # Move back one step to check for new pairs
        elif tokens[i] == 'A#' and tokens[i + 1] == '#B':
            # A# #B becomes #B A#
            tokens[i], tokens[i + 1] = '#B', 'A#'
            i += 1
        elif tokens[i] == 'B#' and tokens[i + 1] == '#A':
            # B# #A becomes #A B#
            tokens[i], tokens[i + 1] = '#A', 'B#'
            i += 1
        elif tokens[i] == 'B#' and tokens[i + 1] == '#B':
            # B# #B becomes nothing
            tokens.pop(i)
            tokens.pop(i)
            i = max(i - 1, 0)  # Move back one step to check for new pairs
        else:
            i += 1
    return tokens

# Initial program
program = ['B#', 'B#', '#A', '#B', 'B#', '#A', '#B', '#B', '#A', 'A#']

# Compute the final state
final_state = compute_program(program)

# Print the final state
print(final_state)